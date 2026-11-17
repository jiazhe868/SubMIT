/* forward_model.c - predicted seismograms, moment tensors and penalized misfit of a multiple-subevent model.

   What it does. sub_init_() reads Par.file, the station lists, the windowed data and the Green's
   functions once. sub_forward_() then turns a set of subevents (centroid time, position, depth,
   duration) into predicted teleseismic P, teleseismic SH and regional three-component seismograms,
   solves the five moment-tensor components of every subevent by damped least squares (with the
   total-moment-tensor constraint), aligns every trace on the data by cross-correlation and returns
   the penalized misfit that finv.f90 samples (README section 4). ffwd.f90 calls the same function
   once after enable_forward_outputs_() to write the waveforms, moment tensors and per-station
   residuals that the figures and the model-based station screen read.

   Why one module. Until v2 the data were read into C structs, copied into ~70 Fortran arguments
   and rebuilt into structs on every forward call, and the forward code existed twice
   (code_SubMIT/sub_forward.c and fwd_SubMIT/sub_forward.c, 1300 lines each). Array sizes that
   differed between the copies once overwrote the SH filter coefficients. The data now
   stay in this file; Fortran passes only the model, and arrays are sized from the actual station
   and window counts instead of compile-time maxima.

   History of measured defects that shaped this file:
   - Teleseismic arrival times were referenced to the depth of subevent E1, but the data are
     aligned on the theoretical arrival from the hypocentre depth (SAC t1, computed by
     add_arrival.sh from the header evdp). When E1's depth moved away from the hypocentre, every
     P arrival shifted by (z1-evdp)*eta_p and every SH arrival by (z1-evdp)*eta_s, out of step with
     the absolute timing of the regional Green's functions (cohort report 2026-10-23, item 1: E1
     at 79.8 km for a 110-km hypocentre moved all times by 3.2 s). The reference is now the evdp of
     the data headers.
   - The S-leg part (sP) of each teleseismic P Green's function was shifted with the vertical S
     slowness of the DIRECT-S ray parameter, on top of the depth-dependent sP delay that mtel3
     already builds into the Green's function of each depth (measured on the 2024 Chile GFs at
     63 deg: the sP peak moves 0.318 s/km with depth, eta_p+eta_s from the P ray parameter is
     0.324 s/km). A subevent 40 km from the reference depth had its sP about 3 s early. The whole
     P Green's function now gets one shift, and the P and S legs are summed before convolution.
   - Under -Ofast (no NaN semantics) a NaN misfit could be accepted and froze the chain at NaN
     (cohort report item 4: 20 of 96 chains). The misfit is now tested on its IEEE bits.
   - The cross-correlation of an all-zero synthetic (a zero component weight, or a subevent whose
     pulse lies outside the window) divided 0 by 0 and left the lag uninitialized, which then
     indexed memory outside the trace.
   - The regional zero-fill before the origin time used the depth index instead of the component
     index (harmless only because fresh memory is zero).

   Build: make finv / make ffwd (Makefile). Needs SAC's libsac (design, apply) and MPI (node-shared
   Green's functions when SUBMIT_SHM_GF=1 and MPI is running).
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <mpi.h>
#include "sacsubc.h"

/* Butterworth band-pass from SAC's libsac (xapiir): design once, apply zero-phase to every trace */
void design(int order, char *filter_type, char *prototype, double attenuation, double transition_bandwidth,
            double low_corner_hz, double high_corner_hz, double sampling_interval, float *sn, float *sd, int *nsects);
void apply(float *data, int npts, int zero_phase, float *sn, float *sd, int nsects);
/* moment-tensor decomposition (mtdcmp.f): CLVD fraction and both nodal planes in radians */
void mtdcmp_(double *m0, double *mxx, double *mxy, double *mxz, double *myy, double *myz, double *mzz,
             double *clvd_fraction, double strike[2], double dip[2], double rake[2]);

#define degrees_per_radian 57.29578
#define max_subevents 9
#define max_filter_sections 30
#define wave_p 0
#define wave_sh 1
#define wave_regional 2

/* ------------------------------------------------------------------------------------------------
   state read once by sub_init_ and used by every forward call
   ------------------------------------------------------------------------------------------------ */
typedef struct {
    char data_file[256];             /* windowed record, e.g. data/IU.ANMO.00.z or data/vel_IU.ANMO.00.z */
    float azimuth_degrees;           /* event-to-station azimuth */
    float p_apparent_velocity;       /* 111.195 km/deg / TauP P ray parameter (km/s) */
    float s_apparent_velocity;       /* same for S */
    char gf_distance_key[64];        /* distance part of the GF file name, e.g. 8713 or vel_8713 */
} teleseismic_trace;

typedef struct {
    char file_prefix[256];           /* dataloc/NET.STA.LOC. ; the e, n, z records append one letter */
    float longitude, latitude;
    float component_weight[3];       /* e, n, z (0 drops a component cleanly) */
} regional_station;

static struct {
    float sampling_interval;
    float window_begin[3], window_end[3];     /* seconds after t1, per wave type */
    int window_npts[3];
    float stf_begin, stf_end;                 /* moment-rate window, seconds after the origin */
    int stf_npts;
    float filter_sn[3][max_filter_sections], filter_sd[3][max_filter_sections];
    int filter_sections[3];
    float wave_weight[3];
    float max_shift_seconds[3];
    float tikhonov_alpha, cmt_scaling, dc_constraint;
    float total_moment_tensor[5];             /* Mxx Myy Mxy Mxz Myz, 1e27 dyne-cm (totalmt.dat) */
    float epicentre_longitude, epicentre_latitude;
    float source_p_velocity, source_s_velocity;
    float hypocentre_depth_km;                /* evdp of the data headers: depth of the t1 picks */
    float gf_depth_first, gf_depth_step;
    int gf_depth_count;
    float regional_distance_first, regional_distance_step;
    int regional_distance_count;
    int trace_count[3];                       /* P traces (displacement then velocity), SH traces, regional stations */
    teleseismic_trace *p_traces, *sh_traces;
    regional_station *regional_stations;
    float *p_data, *sh_data, *regional_data;  /* filtered windows; regional stores e, n, z per station */
    float *p_green, *sh_green, *regional_green;
    int outputs_enabled;
} setup;

/* Green's-function layout: [depth][trace][component][sample] for P (6 components: dd, ds, ss
   fault types, each split into its P leg and S leg as mtel3 writes them) and SH (dsst, ssst);
   [depth][distance][component][sample] for regional (14: 6 vertical, 6 radial, 2 tangential). */
static const char *p_component_names[6] = {"ddpz", "ddsz", "dspz", "dssz", "sspz", "sssz"};
static const char *sh_component_names[2] = {"dsst", "ssst"};
static const char *regional_component_names[14] = {"ddpz", "ddsz", "dspz", "dssz", "sspz", "sssz",
                                                   "ddpr", "ddsr", "dspr", "dssr", "sspr", "sssr", "dsst", "ssst"};

static float *p_green_at(int depth_index, int trace, int component) {
    return setup.p_green + ((size_t)(depth_index * setup.trace_count[wave_p] + trace) * 6 + component) * setup.window_npts[wave_p];
}
static float *sh_green_at(int depth_index, int trace, int component) {
    return setup.sh_green + ((size_t)(depth_index * setup.trace_count[wave_sh] + trace) * 2 + component) * setup.window_npts[wave_sh];
}
static float *regional_green_at(int depth_index, int distance_index, int component) {
    return setup.regional_green + ((size_t)(depth_index * setup.regional_distance_count + distance_index) * 14 + component)
                                  * setup.window_npts[wave_regional];
}

static void fail(const char *message, const char *detail) {
    fprintf(stderr, "ERROR: %s%s%s\n", message, detail ? " " : "", detail ? detail : "");
    exit(1);
}

/* ------------------------------------------------------------------------------------------------
   reading: Par.file, station lists, SAC windows
   ------------------------------------------------------------------------------------------------ */

/* Keyword parser: "key= value #comment" lines in any order, unknown keys ignored, every missing
   required key fatal and named (a positional reader once misassigned every parameter after one
   deleted line). Optional keys keep their defaults. */
static void read_par_file(char *p_station_file, char *sh_station_file, char *regional_station_file,
                          char *gf_prefix, char *gf_separator, char *gf_suffix, char *total_mt_file) {
    struct par_key { const char *name; float *real_value; int *int_value; char *text_value; int required; int seen; };
    float low_hz[3], high_hz[3], last_depth, last_distance;
    int ignored_count;
    struct par_key keys[] = {
        {"bg_timeP=", &setup.window_begin[wave_p]}, {"nd_timeP=", &setup.window_end[wave_p]},
        {"bg_timeSH=", &setup.window_begin[wave_sh]}, {"nd_timeSH=", &setup.window_end[wave_sh]},
        {"bg_timeRayl=", &setup.window_begin[wave_regional]}, {"nd_timeRayl=", &setup.window_end[wave_regional]},
        {"stf_btime=", &setup.stf_begin}, {"stf_etime=", &setup.stf_end}, {"deltat=", &setup.sampling_interval},
        {"weightP=", &setup.wave_weight[wave_p]}, {"weightSH=", &setup.wave_weight[wave_sh]},
        {"weightRayl=", &setup.wave_weight[wave_regional]},
        {"Low_frequencyP=", &low_hz[wave_p]}, {"High_frequencyP=", &high_hz[wave_p]},
        {"Low_frequencySH=", &low_hz[wave_sh]}, {"High_frequencySH=", &high_hz[wave_sh]},
        {"Low_frequencyRayl=", &low_hz[wave_regional]}, {"High_frequencyRayl=", &high_hz[wave_regional]},
        {"Tikhonov_alpha=", &setup.tikhonov_alpha}, {"CMTscaling=", &setup.cmt_scaling}, {"DCconstrain=", &setup.dc_constraint},
        {"evlo=", &setup.epicentre_longitude}, {"evla=", &setup.epicentre_latitude},
        {"evvp=", &setup.source_p_velocity}, {"evvs=", &setup.source_s_velocity},
        {"bg_localdist=", &setup.regional_distance_first}, {"nd_localdist=", &last_distance},
        {"localdist_interval=", &setup.regional_distance_step},
        {"bg_depth=", &setup.gf_depth_first}, {"nd_depth=", &last_depth}, {"depth_interval=", &setup.gf_depth_step},
        {"num_sta_P=", NULL, &setup.trace_count[wave_p]}, {"num_sta_SH=", NULL, &setup.trace_count[wave_sh]},
        {"num_sta_Rayl=", NULL, &setup.trace_count[wave_regional]},
        {"stainfofileP=", NULL, NULL, p_station_file}, {"stainfofileSH=", NULL, NULL, sh_station_file},
        {"stainfofileRayl=", NULL, NULL, regional_station_file},
        {"greenchar1=", NULL, NULL, gf_prefix}, {"greenchar2=", NULL, NULL, gf_separator}, {"greenchar3=", NULL, NULL, gf_suffix},
        /* optional: legacy runs without these keys keep the historical values */
        {"maxshftP=", &setup.max_shift_seconds[wave_p], NULL, NULL, -1},
        {"maxshftSH=", &setup.max_shift_seconds[wave_sh], NULL, NULL, -1},
        {"maxshftRayl=", &setup.max_shift_seconds[wave_regional], NULL, NULL, -1},
        {"TotalMomentTensor=", NULL, NULL, total_mt_file, -1},
        {"num_subevent=", NULL, &ignored_count, NULL, -1},
    };
    int key_count = sizeof keys / sizeof keys[0], key_index, missing_count = 0, wave;
    char line[1024], name[256], value[512];
    FILE *par_file = fopen("Par.file", "r");
    if (par_file == NULL) fail("cannot open Par.file", NULL);
    setup.max_shift_seconds[wave_p] = 2;       /* historical defaults (s) */
    setup.max_shift_seconds[wave_sh] = 5;
    setup.max_shift_seconds[wave_regional] = 3;
    strcpy(total_mt_file, "totalmt.dat");
    for (key_index = 0; key_index < key_count; key_index++) {
        keys[key_index].required = keys[key_index].required == 0;
        keys[key_index].seen = 0;
    }
    while (fgets(line, sizeof line, par_file)) {
        if (sscanf(line, "%255s %511s", name, value) < 2) continue;
        for (key_index = 0; key_index < key_count; key_index++) {
            struct par_key *key = &keys[key_index];
            if (strcmp(name, key->name) != 0) continue;
            if (key->real_value) sscanf(value, "%f", key->real_value);
            else if (key->int_value) sscanf(value, "%d", key->int_value);
            else { strncpy(key->text_value, value, 255); key->text_value[255] = 0; }
            key->seen = 1;
            break;
        }
    }
    fclose(par_file);
    for (key_index = 0; key_index < key_count; key_index++) {
        if (keys[key_index].required && !keys[key_index].seen) {
            fprintf(stderr, "ERROR: Par.file is missing the required key '%s'\n", keys[key_index].name);
            missing_count++;
        }
    }
    if (missing_count) exit(1);
    for (wave = 0; wave < 3; wave++) {
        char filter_type[] = "BP", prototype[] = "BU";
        setup.window_npts[wave] = (int)((setup.window_end[wave] - setup.window_begin[wave]) / setup.sampling_interval + 1);
        design(4, filter_type, prototype, 1.0, 1.0, (double)low_hz[wave], (double)high_hz[wave],
               (double)setup.sampling_interval, setup.filter_sn[wave], setup.filter_sd[wave], &setup.filter_sections[wave]);
    }
    setup.stf_npts = (int)((setup.stf_end - setup.stf_begin) / setup.sampling_interval + 1);
    if (fabs(fmod((last_depth - setup.gf_depth_first) / setup.gf_depth_step * 100, 100)) > 0.01)
        fail("Par.file: nd_depth - bg_depth is not a multiple of depth_interval", NULL);
    if (fabs(fmod((last_distance - setup.regional_distance_first) / setup.regional_distance_step * 100, 100)) > 0.01)
        fail("Par.file: nd_localdist - bg_localdist is not a multiple of localdist_interval", NULL);
    setup.gf_depth_count = (int)((last_depth - setup.gf_depth_first) / setup.gf_depth_step + 1);
    setup.regional_distance_count = (int)((last_distance - setup.regional_distance_first) / setup.regional_distance_step + 1);
    if (setup.wave_weight[wave_p] <= 0 || setup.wave_weight[wave_sh] <= 0 || setup.wave_weight[wave_regional] <= 0)
        fail("Par.file: weightP, weightSH and weightRayl must all be positive", NULL);
}

static void read_teleseismic_station_list(const char *file_name, int trace_count, teleseismic_trace **traces) {
    float distance_degrees;
    int trace;
    FILE *list = fopen(file_name, "r");
    if (list == NULL) fail("cannot open station list", file_name);
    *traces = calloc(trace_count > 0 ? trace_count : 1, sizeof(teleseismic_trace));
    for (trace = 0; trace < trace_count; trace++) {
        teleseismic_trace *t = &(*traces)[trace];
        if (fscanf(list, "%255s %f %f %f %f %63s", t->data_file, &distance_degrees, &t->azimuth_degrees,
                   &t->p_apparent_velocity, &t->s_apparent_velocity, t->gf_distance_key) != 6) {
            fprintf(stderr, "ERROR: %s has fewer than the %d entries that Par.file announces (stale count?)\n", file_name, trace_count);
            exit(1);
        }
    }
    fclose(list);
}

static void read_regional_station_list(const char *file_name, int station_count) {
    char network_station[64];
    int station;
    FILE *list = fopen(file_name, "r");
    if (list == NULL) fail("cannot open station list", file_name);
    setup.regional_stations = calloc(station_count > 0 ? station_count : 1, sizeof(regional_station));
    for (station = 0; station < station_count; station++) {
        regional_station *s = &setup.regional_stations[station];
        if (fscanf(list, "%255s %f %f %63s %f %f %f", s->file_prefix, &s->longitude, &s->latitude, network_station,
                   &s->component_weight[0], &s->component_weight[1], &s->component_weight[2]) != 7) {
            fprintf(stderr, "ERROR: %s has fewer than the %d entries that Par.file announces (stale count?)\n", file_name, station_count);
            exit(1);
        }
    }
    fclose(list);
}

/* Copy npts samples starting window_begin seconds after the record's t1 (header T1). Samples outside
   the record are zero; the count is returned so callers can refuse windows that miss the record. */
static int read_sac_window(const char *file_name, float window_begin, int npts, float *out, float *hypocentre_depth) {
    float *samples = NULL, begin_time, pick_time, record_interval;
    int error, record_npts, first_index, sample, outside_count = 0;
    brsac(1000000, (char *)file_name, &samples, &error);
    if (error < 0 && error != -2) fail("cannot read SAC file", file_name);
    getfhv("B", &begin_time, &error);
    getfhv("T1", &pick_time, &error);
    getfhv("DELTA", &record_interval, &error);
    getnhv("NPTS", &record_npts, &error);
    if (pick_time == -12345.0f) fail("SAC header t1 (arrival pick) is not set in", file_name);
    if (fabsf(record_interval - setup.sampling_interval) > 1e-4f * setup.sampling_interval)
        fail("sampling interval differs from Par.file deltat in", file_name);
    if (hypocentre_depth) getfhv("EVDP", hypocentre_depth, &error);
    first_index = (int)((pick_time + window_begin - begin_time) / setup.sampling_interval);
    for (sample = 0; sample < npts; sample++) {
        int index = first_index + sample;
        if (index >= 0 && index < record_npts) out[sample] = samples[index];
        else { out[sample] = 0; outside_count++; }
    }
    free(samples);
    return outside_count;
}

static void read_data(void) {
    int wave, trace, component, npts;
    char file_name[300];
    const char *component_letters = "enz";
    setup.hypocentre_depth_km = -1;
    for (wave = wave_p; wave <= wave_sh; wave++) {
        teleseismic_trace *traces = wave == wave_p ? setup.p_traces : setup.sh_traces;
        int count = setup.trace_count[wave];
        float *data;
        npts = setup.window_npts[wave];
        data = calloc((size_t)(count > 0 ? count : 1) * npts, sizeof(float));
        for (trace = 0; trace < count; trace++) {
            float depth = -1;
            float *window = data + (size_t)trace * npts;
            if (read_sac_window(traces[trace].data_file, setup.window_begin[wave], npts, window, &depth) > npts / 2)
                fail("the inversion window lies mostly outside the record", traces[trace].data_file);
            if (setup.hypocentre_depth_km < 0 && depth >= 0 && depth != -12345.0f) setup.hypocentre_depth_km = depth;
            apply(window, npts, 1, setup.filter_sn[wave], setup.filter_sd[wave], setup.filter_sections[wave]);
            /* the second half of the P list are velocity records, weighted x2 (their GFs too) */
            if (wave == wave_p && trace >= count / 2) {
                int sample;
                for (sample = 0; sample < npts; sample++) window[sample] *= 2;
            }
        }
        if (wave == wave_p) setup.p_data = data; else setup.sh_data = data;
    }
    npts = setup.window_npts[wave_regional];
    setup.regional_data = calloc((size_t)(setup.trace_count[wave_regional] > 0 ? setup.trace_count[wave_regional] : 1) * 3 * npts, sizeof(float));
    for (trace = 0; trace < setup.trace_count[wave_regional]; trace++) {
        regional_station *station = &setup.regional_stations[trace];
        for (component = 0; component < 3; component++) {
            float *window = setup.regional_data + ((size_t)trace * 3 + component) * npts;
            int sample;
            snprintf(file_name, sizeof file_name, "%s%c", station->file_prefix, component_letters[component]);
            if (read_sac_window(file_name, setup.window_begin[wave_regional], npts, window, NULL) > npts / 2)
                fail("the inversion window lies mostly outside the record", file_name);
            apply(window, npts, 1, setup.filter_sn[wave_regional], setup.filter_sd[wave_regional], setup.filter_sections[wave_regional]);
            for (sample = 0; sample < npts; sample++) window[sample] *= station->component_weight[component];
        }
    }
    if (setup.hypocentre_depth_km < 0) {
        fprintf(stderr, "WARNING: no evdp in the teleseismic data headers - timing referenced to the depth of E1\n");
    }
}

/* ------------------------------------------------------------------------------------------------
   Green's functions: file names, reading, optional node-shared storage
   ------------------------------------------------------------------------------------------------ */

/* depth or distance as it appears in the GF folder/file names: "20" or "20.5" */
static void format_grid_value(float value, char *out) {
    if (fmod(value, 1) < 0.01) sprintf(out, "%d", (int)value);
    else if (fabs(fmod(value * 100 + 0.01, 10) / 100) < 0.01) sprintf(out, "%.1f", value);
    else fail("GF depths and distances may carry at most one decimal; change the grid in Par.file", NULL);
}

static void read_green_file(const char *path, float window_begin, int npts, float *out, float scale) {
    int sample;
    if (access(path, R_OK) != 0) {
        /* a missing file once left the kernel as silent garbage for one station (2024 Chile AU.MAW) */
        fprintf(stderr, "ERROR: missing Green's function file: %s\n", path);
        fprintf(stderr, "       regenerate the Green's functions or remove the station from the station lists\n");
        exit(1);
    }
    read_sac_window(path, window_begin, npts, out, NULL);
    for (sample = 0; sample < npts; sample++) out[sample] *= scale;
}

static void read_green_functions(const char *prefix, const char *separator, const char *suffix) {
    int depth_index, trace, component, distance_index;
    char depth_text[32], distance_text[32], path[1024];
    for (depth_index = 0; depth_index < setup.gf_depth_count; depth_index++) {
        format_grid_value(setup.gf_depth_first + depth_index * setup.gf_depth_step, depth_text);
        for (trace = 0; trace < setup.trace_count[wave_p]; trace++) {
            /* velocity traces (second half of the list) are weighted x2 like their data */
            float scale = trace < setup.trace_count[wave_p] / 2 ? 1e5f : 2e5f;
            for (component = 0; component < 6; component++) {
                snprintf(path, sizeof path, "%s%s%s%s%s%s", prefix, depth_text, separator,
                         setup.p_traces[trace].gf_distance_key, suffix, p_component_names[component]);
                read_green_file(path, setup.window_begin[wave_p], setup.window_npts[wave_p], p_green_at(depth_index, trace, component), scale);
            }
        }
        for (trace = 0; trace < setup.trace_count[wave_sh]; trace++) {
            for (component = 0; component < 2; component++) {
                snprintf(path, sizeof path, "%s%s%s%s%s%s", prefix, depth_text, separator,
                         setup.sh_traces[trace].gf_distance_key, suffix, sh_component_names[component]);
                read_green_file(path, setup.window_begin[wave_sh], setup.window_npts[wave_sh], sh_green_at(depth_index, trace, component), 1e5f);
            }
        }
        for (distance_index = 0; distance_index < setup.regional_distance_count; distance_index++) {
            int npts = setup.window_npts[wave_regional];
            /* fk records may start after bg_timeRayl: read from the origin time (t1 = 0) and keep the
               samples before it at zero */
            int samples_before_origin = (int)((0 - setup.window_begin[wave_regional]) / setup.sampling_interval);
            if (samples_before_origin < 0 || samples_before_origin >= npts) fail("bg_timeRayl must lie before the origin time", NULL);
            format_grid_value(setup.regional_distance_first + distance_index * setup.regional_distance_step, distance_text);
            for (component = 0; component < 14; component++) {
                float *green = regional_green_at(depth_index, distance_index, component);
                snprintf(path, sizeof path, "%s%s%s%s%s%s", prefix, depth_text, separator, distance_text, suffix,
                         regional_component_names[component]);
                memset(green, 0, sizeof(float) * samples_before_origin);
                read_green_file(path, 0, npts - samples_before_origin, green + samples_before_origin, 1e5f);
            }
        }
    }
}

/* SUBMIT_SHM_GF=1 under MPI: one physical copy of the Green's functions per node in POSIX shared
   memory, read from disk once (validated bit-identical to per-rank copies; 39% faster at 8 ranks on
   Calama and ~20 GB of RAM saved at 32 ranks). Falls back to per-rank memory when /dev/shm is too
   small (posix_fallocate fails early instead of SIGBUS mid-fill). */
static void allocate_and_read_green_functions(const char *prefix, const char *separator, const char *suffix) {
    size_t p_floats = (size_t)setup.gf_depth_count * setup.trace_count[wave_p] * 6 * setup.window_npts[wave_p];
    size_t sh_floats = (size_t)setup.gf_depth_count * setup.trace_count[wave_sh] * 2 * setup.window_npts[wave_sh];
    size_t regional_floats = (size_t)setup.gf_depth_count * setup.regional_distance_count * 14 * setup.window_npts[wave_regional];
    size_t total_bytes = (p_floats + sh_floats + regional_floats) * sizeof(float);
    const char *shared_request = getenv("SUBMIT_SHM_GF");
    int mpi_running = 0;
    if (shared_request && shared_request[0] == '1') MPI_Initialized(&mpi_running);
    if (mpi_running) {
        MPI_Comm node_communicator;
        int node_rank, descriptor = -1;
        long segment_tag = 0;
        char segment_name[128];
        void *base;
        MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &node_communicator);
        MPI_Comm_rank(node_communicator, &node_rank);
        if (node_rank == 0) segment_tag = (long)getpid();
        MPI_Bcast(&segment_tag, 1, MPI_LONG, 0, node_communicator);
        snprintf(segment_name, sizeof segment_name, "/submit_gf_%ld", segment_tag);
        if (node_rank == 0) {
            descriptor = shm_open(segment_name, O_CREAT | O_RDWR, 0600);
            if (descriptor >= 0 && (ftruncate(descriptor, (off_t)total_bytes) != 0 ||
                                    posix_fallocate(descriptor, 0, (off_t)total_bytes) != 0)) {
                close(descriptor); shm_unlink(segment_name); descriptor = -1;
                fprintf(stderr, "SUBMIT_SHM_GF: /dev/shm too small for %.1f GB - per-rank copies\n", total_bytes / 1073741824.0);
            }
        }
        MPI_Barrier(node_communicator);
        if (node_rank != 0) descriptor = shm_open(segment_name, O_RDWR, 0600);
        base = descriptor >= 0 ? mmap(NULL, total_bytes, PROT_READ | PROT_WRITE, MAP_SHARED, descriptor, 0) : MAP_FAILED;
        if (descriptor >= 0) close(descriptor);
        MPI_Barrier(node_communicator);
        if (node_rank == 0) shm_unlink(segment_name);   /* the segment vanishes when the last rank exits */
        /* every rank must agree, otherwise some would wait for a fill that never comes */
        {
            int mapped_here = base != MAP_FAILED, mapped_everywhere = 0;
            MPI_Allreduce(&mapped_here, &mapped_everywhere, 1, MPI_INT, MPI_MIN, node_communicator);
            if (mapped_everywhere) {
                setup.p_green = base;
                setup.sh_green = setup.p_green + p_floats;
                setup.regional_green = setup.sh_green + sh_floats;
                if (node_rank == 0) {
                    read_green_functions(prefix, separator, suffix);
                    fprintf(stderr, "SUBMIT_SHM_GF: node-shared Green's functions, one %.1f GB copy\n", total_bytes / 1073741824.0);
                }
                MPI_Barrier(node_communicator);
                MPI_Comm_free(&node_communicator);
                return;
            }
            if (base != MAP_FAILED) munmap(base, total_bytes);
            if (node_rank == 0) fprintf(stderr, "SUBMIT_SHM_GF: shared memory unavailable - per-rank copies\n");
        }
        MPI_Comm_free(&node_communicator);
    }
    setup.p_green = malloc((p_floats + sh_floats + regional_floats + 1) * sizeof(float));
    if (setup.p_green == NULL) fail("not enough memory for the Green's functions", NULL);
    setup.sh_green = setup.p_green + p_floats;
    setup.regional_green = setup.sh_green + sh_floats;
    read_green_functions(prefix, separator, suffix);
}

/* ------------------------------------------------------------------------------------------------
   Fortran entry points: initialisation
   ------------------------------------------------------------------------------------------------ */

/* Read everything once. Returns the S velocity at the source (finv's causality cap is a multiple). */
void sub_init_(float *source_s_velocity) {
    char p_list[256], sh_list[256], regional_list[256], gf_prefix[256], gf_separator[256], gf_suffix[256], total_mt_file[256];
    FILE *total_mt;
    int component;
    read_par_file(p_list, sh_list, regional_list, gf_prefix, gf_separator, gf_suffix, total_mt_file);
    read_teleseismic_station_list(p_list, setup.trace_count[wave_p], &setup.p_traces);
    read_teleseismic_station_list(sh_list, setup.trace_count[wave_sh], &setup.sh_traces);
    read_regional_station_list(regional_list, setup.trace_count[wave_regional]);
    total_mt = fopen(total_mt_file, "r");
    if (total_mt == NULL) fail("cannot open the total moment tensor", total_mt_file);
    for (component = 0; component < 5; component++)
        if (fscanf(total_mt, "%f", &setup.total_moment_tensor[component]) != 1) fail("totalmt.dat needs 5 numbers", NULL);
    fclose(total_mt);
    read_data();
    allocate_and_read_green_functions(gf_prefix, gf_separator, gf_suffix);
    *source_s_velocity = setup.source_s_velocity;
}

void enable_forward_outputs_(void) {
    setup.outputs_enabled = 1;
}

/* ------------------------------------------------------------------------------------------------
   synthetics: moment-rate pulse, convolution, basis seismograms for the 5 moment-tensor components
   ------------------------------------------------------------------------------------------------ */

/* Gaussian moment-rate pulse of unit area whose full width at 10% of the peak is apparent_duration,
   centred at centre_time, truncated at +/- apparent_duration/1.5 like the published runs. Fills
   rate[first..last-1] on the moment-rate grid stf_begin + m*dt and returns first/last. */
static void gaussian_moment_rate(float centre_time, float apparent_duration, float *rate, int *first, int *last) {
    float dt = setup.sampling_interval;
    float sigma = apparent_duration / (2 * (float)sqrt(2 * log(10)));
    float variance_twice = sigma * sigma * 2.0f;
    int sample;
    *first = (int)((centre_time - apparent_duration / 1.5f - setup.stf_begin) / dt);
    *last = (int)((centre_time + apparent_duration / 1.5f - setup.stf_begin) / dt);
    if (*first < 0) *first = 0;
    if (*first > setup.stf_npts) *first = setup.stf_npts;
    if (*last > setup.stf_npts - 1) *last = setup.stf_npts - 1;
    if (*last < 0) *last = 0;
    for (sample = *first; sample < *last; sample++) {
        float time = setup.stf_begin + sample * dt;
        rate[sample] = dt / ((float)sqrt(2 * 3.1416) * sigma) * expf(-(time - centre_time) * (time - centre_time) / variance_twice);
    }
}

/* out[q] = sum_m rate[m] green[q + offset - m] over the nonzero part of the pulse only: the window of
   the data (npts samples starting window_begin after t1) shifted by the moment-rate window start. */
static void convolve_with_moment_rate(const float *rate, int first, int last, const float *green, int npts, int offset, float *out) {
    int q, m;
    for (q = 0; q < npts; q++) {
        int m_low = q + offset - npts + 1 > first ? q + offset - npts + 1 : first;
        int m_high = q + offset + 1 < last ? q + offset + 1 : last;
        float sum = 0;
        for (m = m_low; m < m_high; m++) sum += rate[m] * green[q + offset - m];
        out[q] = sum;
    }
}

static int convolution_offset(int wave) {
    /* index of the first data-window sample in the full convolution (as the published code computed it) */
    float data_begin = setup.window_begin[wave];
    return (int)((data_begin - (data_begin + setup.stf_begin)) / setup.sampling_interval + 1) - 1;
}

/* the five basis seismograms of one vertical-type (P or regional Z) record from the dd, ds, ss parts */
static void vertical_basis(const float *dd, const float *ds, const float *ss, float azimuth_radians, int npts, float *basis[5]) {
    float cos2 = cosf(2 * azimuth_radians), sin2 = sinf(2 * azimuth_radians), cos1 = cosf(azimuth_radians), sin1 = sinf(azimuth_radians);
    int q;
    for (q = 0; q < npts; q++) {
        basis[0][q] = -ss[q] / 2 * cos2 - dd[q] / 2;
        basis[1][q] = ss[q] / 2 * cos2 - dd[q] / 2;
        basis[2][q] = -ss[q] * sin2;
        basis[3][q] = -ds[q] * cos1;
        basis[4][q] = -ds[q] * sin1;
    }
}

/* the five basis seismograms of one transverse (SH or regional T) record from the ds, ss parts */
static void transverse_basis(const float *ds, const float *ss, float azimuth_radians, int npts, float *basis[5]) {
    float cos2 = cosf(2 * azimuth_radians), sin2 = sinf(2 * azimuth_radians), cos1 = cosf(azimuth_radians), sin1 = sinf(azimuth_radians);
    int q;
    for (q = 0; q < npts; q++) {
        basis[0][q] = -ss[q] / 2 * sin2;
        basis[1][q] = ss[q] / 2 * sin2;
        basis[2][q] = ss[q] * cos2;
        basis[3][q] = -ds[q] * sin1;
        basis[4][q] = ds[q] * cos1;
    }
}

/* ------------------------------------------------------------------------------------------------
   alignment and misfit
   ------------------------------------------------------------------------------------------------ */

/* Best normalized cross-correlation of the synthetic against the data over lags -max_lag..max_lag-1
   samples (the search of the published code). A positive lag means the synthetic arrives later:
   synthetic[k] matches data[k - lag]. An all-zero trace returns lag 0 and correlation 0 instead of
   0/0. */
static void align_by_cross_correlation(const float *data, const float *synthetic, int npts, int max_lag, int *best_lag, float *best_correlation) {
    int overlap = npts - max_lag, shift_index, q;
    double data_energy = 0, synthetic_total = 0;
    *best_lag = 0;
    *best_correlation = 0;
    if (overlap <= 0) return;
    for (q = 0; q < overlap; q++) data_energy += (double)data[q + max_lag] * data[q + max_lag];
    for (q = 0; q < npts; q++) synthetic_total += (double)synthetic[q] * synthetic[q];
    if (data_energy <= 0 || synthetic_total <= 0) return;
    *best_correlation = -2;
    for (shift_index = 0; shift_index < 2 * max_lag; shift_index++) {
        double cross = 0, synthetic_energy = 0, correlation;
        for (q = 0; q < overlap; q++) {
            int synthetic_index = q + shift_index;
            float s = synthetic_index < npts ? synthetic[synthetic_index] : 0;
            cross += (double)data[q + max_lag] * s;
            synthetic_energy += (double)s * s;
        }
        correlation = synthetic_energy > 0 ? cross / sqrt(data_energy * synthetic_energy) : 0;
        if (correlation > *best_correlation) {
            *best_correlation = (float)correlation;
            *best_lag = shift_index - max_lag;
        }
    }
    if (max_lag == 0) *best_correlation = 0;
}

/* NaN/Inf test on the IEEE-754 bit pattern: -Ofast implies -ffinite-math-only, under which the
   compiler may fold isfinite() to "true" (observed with Intel icx 2025.3, cohort report item 4). */
static int is_finite_float(float value) {
    unsigned int bits;
    memcpy(&bits, &value, sizeof bits);
    return ((bits >> 23) & 0xffu) != 0xffu;
}

/* Cholesky factorization in place of a symmetric positive-definite n x n matrix (row-major, lower
   triangle used); returns 0 when a pivot is not positive. */
static int cholesky_factorize(double *matrix, int n) {
    int i, j, k;
    for (j = 0; j < n; j++) {
        double pivot = matrix[j * n + j];
        for (k = 0; k < j; k++) pivot -= matrix[j * n + k] * matrix[j * n + k];
        if (!(pivot > 0)) return 0;
        matrix[j * n + j] = sqrt(pivot);
        for (i = j + 1; i < n; i++) {
            double value = matrix[i * n + j];
            for (k = 0; k < j; k++) value -= matrix[i * n + k] * matrix[j * n + k];
            matrix[i * n + j] = value / matrix[j * n + j];
        }
    }
    return 1;
}

static void cholesky_solve(const double *factor, int n, const double *rhs, double *solution) {
    int i, k;
    for (i = 0; i < n; i++) {
        double value = rhs[i];
        for (k = 0; k < i; k++) value -= factor[i * n + k] * solution[k];
        solution[i] = value / factor[i * n + i];
    }
    for (i = n - 1; i >= 0; i--) {
        double value = solution[i];
        for (k = i + 1; k < n; k++) value -= factor[k * n + i] * solution[k];
        solution[i] = value / factor[i * n + i];
    }
}

/* P and T axes of a double couple (strike, dip, rake in radians; x north, y east, z down). The two
   labels come out swapped relative to Aki & Richards, which does not matter: the stress penalty
   treats both axes alike. */
static void pressure_and_tension_axes(double strike, double dip, double rake, double p_axis[3], double t_axis[3]) {
    double normal[3] = {-sin(dip) * sin(strike), sin(dip) * cos(strike), -cos(dip)};
    double slip[3] = {cos(rake) * cos(strike) + sin(rake) * cos(dip) * sin(strike),
                      cos(rake) * sin(strike) - sin(rake) * cos(dip) * cos(strike), -sin(rake) * sin(dip)};
    double p_norm = 0, t_norm = 0;
    int i;
    for (i = 0; i < 3; i++) {
        p_axis[i] = (normal[i] + slip[i]) / 2;
        t_axis[i] = (normal[i] - slip[i]) / 2;
        p_norm += p_axis[i] * p_axis[i];
        t_norm += t_axis[i] * t_axis[i];
    }
    p_norm = sqrt(p_norm) > 0 ? sqrt(p_norm) : 1e-6;
    t_norm = sqrt(t_norm) > 0 ? sqrt(t_norm) : 1e-6;
    for (i = 0; i < 3; i++) { p_axis[i] /= p_norm; t_axis[i] /= t_norm; }
}

static double angle_between_axes_degrees(const double first[3], const double second[3]) {
    double cosine = first[0] * second[0] + first[1] * second[1] + first[2] * second[2];
    double angle;
    if (cosine >= 1) cosine -= 1e-6;
    if (cosine <= -1) cosine += 1e-6;
    angle = acos(cosine) * degrees_per_radian;
    return angle > 90 ? 180 - angle : angle;      /* axes have no sign */
}

/* ------------------------------------------------------------------------------------------------
   outputs of ffwd (figures, station screen, weight balancing)
   ------------------------------------------------------------------------------------------------ */

static void write_sac_trace(const char *file_name, const float *samples, int npts) {
    float maximum, minimum, mean;
    int error;
    scmxmn((float *)samples, npts, &maximum, &minimum, &mean);
    newhdr();
    setfhv("DEPMAX", maximum, &error);
    setfhv("DEPMIN", minimum, &error);
    setfhv("DEPMEN", mean, &error);
    setnhv("NPTS    ", npts, &error);
    setfhv("DELTA   ", setup.sampling_interval, &error);
    setfhv("B       ", 0, &error);
    setihv("IFTYPE  ", "ITIME   ", &error);
    setfhv("E       ", (npts - 1) * setup.sampling_interval, &error);
    setlhv("LEVEN   ", 1, &error);
    setlhv("LOVROK  ", 1, &error);
    setlhv("LCALDA  ", 1, &error);
    setnhv("NZYEAR", 1970, &error);
    setnhv("NZJDAY", 1, &error);
    setnhv("NZHOUR", 0, &error);
    setnhv("NZMIN", 0, &error);
    setnhv("NZSEC", 0, &error);
    setnhv("NZMSEC", 0, &error);
    bwsac(npts, (char *)file_name, (float *)samples);
}

static void remove_sac_files(const char *folder) {
    char command[512];
    snprintf(command, sizeof command, "mkdir -p %s && rm -f %s/*.sac", folder, folder);
    if (system(command) != 0) fprintf(stderr, "WARNING: could not clean %s\n", folder);
}

/* ------------------------------------------------------------------------------------------------
   Fortran entry point: one forward evaluation
   ------------------------------------------------------------------------------------------------ */

/* model: 7 numbers per subevent, in the order of Input.model and finv's convertModel:
   centroid time (s), x east (km), y north (km), duration (s), rupture velocity (km/s),
   rupture direction (deg), centroid depth (km). Returns the penalized misfit, 9999 when it is not
   finite or when no moment could be solved. */
void sub_forward_(int *subevent_count_pointer, float *model, float *penalized_misfit) {
    const int subevent_count = *subevent_count_pointer;
    const int unknown_count = 5 * subevent_count;
    const float dt = setup.sampling_interval;
    const int npts_p = setup.window_npts[wave_p], npts_sh = setup.window_npts[wave_sh], npts_regional = setup.window_npts[wave_regional];
    const int p_count = setup.trace_count[wave_p], sh_count = setup.trace_count[wave_sh], regional_count = setup.trace_count[wave_regional];
    /* rows of the data vector: P traces, then SH traces, then e, n, z of every regional station */
    const int trace_total = p_count + sh_count + 3 * regional_count;
    const size_t row_count = (size_t)p_count * npts_p + (size_t)sh_count * npts_sh + (size_t)3 * regional_count * npts_regional;
    const float reference_depth = setup.hypocentre_depth_km >= 0 ? setup.hypocentre_depth_km : model[6];
    int *trace_npts = malloc(sizeof(int) * (trace_total > 0 ? trace_total : 1));
    size_t *trace_start = malloc(sizeof(size_t) * (trace_total > 0 ? trace_total : 1));
    int *trace_wave = malloc(sizeof(int) * (trace_total > 0 ? trace_total : 1));
    float *data_vector = malloc(sizeof(float) * (row_count + 1));
    float *shifted_data = malloc(sizeof(float) * (row_count + 1));
    float *synthetic = malloc(sizeof(float) * (row_count + 1));
    float *design_matrix = calloc((size_t)unknown_count * row_count + 1, sizeof(float));   /* column per unknown */
    float *moment_rate = calloc(setup.stf_npts + 1, sizeof(float));
    int max_npts = npts_p > npts_sh ? npts_p : npts_sh;
    float *work = NULL, *basis[5], *regional_basis[3][5];
    double normal_matrix[45 * 45], rhs[45], solution[45], moment_tensor[45];
    double data_energy = 0, cmt_weight, moment_cmt;
    double residual_sum[3] = {0, 0, 0}, norm_sum[3] = {0, 0, 0};
    float *trace_correlation = malloc(sizeof(float) * (trace_total > 0 ? trace_total : 1));
    float *trace_l2 = malloc(sizeof(float) * (trace_total > 0 ? trace_total : 1));
    int *trace_lag = malloc(sizeof(int) * (trace_total > 0 ? trace_total : 1));
    int subevent, trace, q, column, pass, i, k;
    size_t row;
    if (npts_regional > max_npts) max_npts = npts_regional;
    work = malloc(sizeof(float) * (size_t)max_npts * 30);
    for (k = 0; k < 5; k++) basis[k] = work + (size_t)k * max_npts;
    for (i = 0; i < 3; i++) for (k = 0; k < 5; k++) regional_basis[i][k] = work + (size_t)(5 + 5 * i + k) * max_npts;
    *penalized_misfit = 9999.0f;
    if (subevent_count < 1 || subevent_count > max_subevents) fail("subevent count outside 1..9", NULL);

    /* bookkeeping of the trace order and the weighted data vector */
    row = 0;
    for (trace = 0; trace < trace_total; trace++) {
        int wave = trace < p_count ? wave_p : (trace < p_count + sh_count ? wave_sh : wave_regional);
        const float *source;
        trace_wave[trace] = wave;
        trace_npts[trace] = setup.window_npts[wave];
        trace_start[trace] = row;
        if (wave == wave_p) source = setup.p_data + (size_t)trace * npts_p;
        else if (wave == wave_sh) source = setup.sh_data + (size_t)(trace - p_count) * npts_sh;
        else source = setup.regional_data + (size_t)(trace - p_count - sh_count) * npts_regional;
        for (q = 0; q < trace_npts[trace]; q++) {
            data_vector[row + q] = source[q] * setup.wave_weight[wave];
            data_energy += (double)data_vector[row + q] * data_vector[row + q];
        }
        row += trace_npts[trace];
    }

    /* --- columns of the design matrix: one basis seismogram per subevent and moment-tensor component */
    for (subevent = 0; subevent < subevent_count; subevent++) {
        const float *parameters = model + 7 * subevent;
        const float centroid_time = parameters[0], east = parameters[1], north = parameters[2];
        const float duration = parameters[3], rupture_velocity = parameters[4], rupture_direction = parameters[5] / degrees_per_radian;
        const float depth = parameters[6];
        const float offset_km = sqrtf(east * east + north * north);
        const float offset_azimuth = (float)atan2((double)east, (double)north);
        /* linear interpolation between the two grid depths around the subevent; outside the grid
           the nearest grid depth is used */
        int depth_first = (int)floorf((depth - setup.gf_depth_first) / setup.gf_depth_step), depth_second;
        float depth_fraction;
        if (depth_first < 0) depth_first = 0;
        if (depth_first > setup.gf_depth_count - 1) depth_first = setup.gf_depth_count - 1;
        depth_second = depth_first + 1 < setup.gf_depth_count ? depth_first + 1 : depth_first;
        depth_fraction = (depth - (setup.gf_depth_first + depth_first * setup.gf_depth_step)) / setup.gf_depth_step;
        if (depth_second == depth_first || depth_fraction < 0) depth_fraction = 0;
        if (depth_fraction > 1) depth_fraction = 1;

        /* teleseismic P and SH: Green's functions at the station's distance from the epicentre, the
           subevent's position entering as a plane-wave time shift (README eq. for delta t_ij) */
        for (trace = 0; trace < p_count + sh_count; trace++) {
            const int is_p = trace < p_count;
            const teleseismic_trace *station = is_p ? &setup.p_traces[trace] : &setup.sh_traces[trace - p_count];
            const int local = is_p ? trace : trace - p_count, npts = is_p ? npts_p : npts_sh, wave = is_p ? wave_p : wave_sh;
            const float apparent_velocity = is_p ? station->p_apparent_velocity : station->s_apparent_velocity;
            const float source_velocity = is_p ? setup.source_p_velocity : setup.source_s_velocity;
            const float azimuth = station->azimuth_degrees / degrees_per_radian;
            /* vertical slowness at the source with the ray parameter of THIS phase; the S-leg (sP) of
               the P Green's functions needs no separate shift: mtel3 already delays it correctly for
               each grid depth, so the whole P Green's function moves with eta_p */
            const float vertical_slowness = (float)sqrt(1 / (source_velocity * source_velocity) - 1 / (apparent_velocity * apparent_velocity));
            const float apparent_duration = duration - duration * rupture_velocity / apparent_velocity * (float)cos(rupture_direction - azimuth);
            const float centre = centroid_time - offset_km / fabsf(apparent_velocity) * (float)cos(azimuth - offset_azimuth)
                                 - (depth - reference_depth) * vertical_slowness * (apparent_velocity > 0 ? 1 : -1);
            const int offset = convolution_offset(wave);
            float *summed = work + (size_t)20 * max_npts, *convolved = work + (size_t)23 * max_npts;
            int first, last, part;
            gaussian_moment_rate(centre, apparent_duration, moment_rate, &first, &last);
            if (is_p) {
                /* dd, ds, ss: P leg + S leg, interpolated between the two grid depths */
                for (part = 0; part < 3; part++) {
                    const float *a0 = p_green_at(depth_first, local, 2 * part), *b0 = p_green_at(depth_first, local, 2 * part + 1);
                    const float *a1 = p_green_at(depth_second, local, 2 * part), *b1 = p_green_at(depth_second, local, 2 * part + 1);
                    float *sum = summed + (size_t)part * max_npts;
                    for (q = 0; q < npts; q++) sum[q] = (a0[q] + b0[q]) * (1 - depth_fraction) + (a1[q] + b1[q]) * depth_fraction;
                    convolve_with_moment_rate(moment_rate, first, last, sum, npts, offset, convolved + (size_t)part * max_npts);
                }
                vertical_basis(convolved, convolved + max_npts, convolved + 2 * (size_t)max_npts, azimuth, npts, basis);
            } else {
                for (part = 0; part < 2; part++) {
                    const float *a0 = sh_green_at(depth_first, local, part), *a1 = sh_green_at(depth_second, local, part);
                    float *sum = summed + (size_t)part * max_npts;
                    for (q = 0; q < npts; q++) sum[q] = a0[q] * (1 - depth_fraction) + a1[q] * depth_fraction;
                    convolve_with_moment_rate(moment_rate, first, last, sum, npts, offset, convolved + (size_t)part * max_npts);
                }
                transverse_basis(convolved, convolved + max_npts, azimuth, npts, basis);
            }
            for (k = 0; k < 5; k++) {
                float *target = design_matrix + (size_t)(5 * subevent + k) * row_count + trace_start[trace];
                apply(basis[k], npts, 1, setup.filter_sn[wave], setup.filter_sd[wave], setup.filter_sections[wave]);
                for (q = 0; q < npts; q++) target[q] = basis[k][q] * setup.wave_weight[wave];
            }
        }

        /* regional: each subevent uses the Green's functions of its own distance and depth (bilinear),
           so the timing is exact and no shift is needed */
        for (trace = 0; trace < regional_count; trace++) {
            const regional_station *station = &setup.regional_stations[trace];
            const int npts = npts_regional, offset = convolution_offset(wave_regional);
            /* great-circle distance and azimuth from the SUBEVENT (not the epicentre) to the station */
            const float source_longitude = setup.epicentre_longitude + east / (111.32f * (float)cos(setup.epicentre_latitude / degrees_per_radian));
            const float source_latitude = setup.epicentre_latitude + north / 110.574f;
            const double latitude_difference = (station->latitude - source_latitude) / degrees_per_radian;
            const double longitude_difference = (station->longitude - source_longitude) / degrees_per_radian;
            const double haversine = pow(sin(latitude_difference / 2), 2)
                                     + cos(station->latitude / degrees_per_radian) * cos(source_latitude / degrees_per_radian) * pow(sin(longitude_difference / 2), 2);
            const float distance_km = (float)(6371 * 2 * atan2(sqrt(haversine), sqrt(1 - haversine)));
            const float azimuth_radians = (float)atan2(sin(longitude_difference) * cos(station->latitude / degrees_per_radian),
                                                       cos(source_latitude / degrees_per_radian) * sin(station->latitude / degrees_per_radian)
                                                       - sin(source_latitude / degrees_per_radian) * cos(station->latitude / degrees_per_radian) * cos(longitude_difference));
            /* directivity of a finite subevent seen at 0.85 Vs, the regional apparent velocity */
            const float apparent_duration = duration - duration * rupture_velocity / (setup.source_s_velocity * 0.85f) * (float)cos(rupture_direction - azimuth_radians);
            int distance_first = (int)floorf((distance_km - setup.regional_distance_first) / setup.regional_distance_step), distance_second;
            float distance_fraction, node_weight[4];
            float *summed = work + (size_t)20 * max_npts, *convolved = work + (size_t)22 * max_npts;
            int first, last, part, node;
            int component_index[8][2] = {{0, 1}, {2, 3}, {4, 5}, {6, 7}, {8, 9}, {10, 11}, {12, -1}, {13, -1}};  /* Zdd Zds Zss Rdd Rds Rss Tds Tss */
            float *combined[8];
            if (distance_first < 0) distance_first = 0;
            if (distance_first > setup.regional_distance_count - 1) distance_first = setup.regional_distance_count - 1;
            distance_second = distance_first + 1 < setup.regional_distance_count ? distance_first + 1 : distance_first;
            distance_fraction = distance_second == distance_first ? 0 :
                (distance_km - (setup.regional_distance_first + distance_first * setup.regional_distance_step)) / setup.regional_distance_step;
            if (distance_fraction < 0) distance_fraction = 0;
            if (distance_fraction > 1) distance_fraction = 1;
            node_weight[0] = (1 - distance_fraction) * (1 - depth_fraction);
            node_weight[1] = distance_fraction * (1 - depth_fraction);
            node_weight[2] = (1 - distance_fraction) * depth_fraction;
            node_weight[3] = distance_fraction * depth_fraction;
            gaussian_moment_rate(centroid_time, apparent_duration, moment_rate, &first, &last);
            for (part = 0; part < 8; part++) {
                float *sum = summed;
                const int node_depth[4] = {depth_first, depth_first, depth_second, depth_second};
                const int node_distance[4] = {distance_first, distance_second, distance_first, distance_second};
                for (q = 0; q < npts; q++) sum[q] = 0;
                for (node = 0; node < 4; node++) {
                    const float *a = regional_green_at(node_depth[node], node_distance[node], component_index[part][0]);
                    const float *b = component_index[part][1] >= 0 ? regional_green_at(node_depth[node], node_distance[node], component_index[part][1]) : NULL;
                    for (q = 0; q < npts; q++) sum[q] += (a[q] + (b ? b[q] : 0)) * node_weight[node];
                }
                combined[part] = convolved + (size_t)part * max_npts;
                convolve_with_moment_rate(moment_rate, first, last, sum, npts, offset, combined[part]);
            }
            {
                float *radial[5], *transverse[5];
                float sin_az = sinf(azimuth_radians), cos_az = cosf(azimuth_radians);
                int c;
                for (k = 0; k < 5; k++) { radial[k] = regional_basis[0][k]; transverse[k] = regional_basis[1][k]; }
                vertical_basis(combined[3], combined[4], combined[5], azimuth_radians, npts, radial);
                transverse_basis(combined[6], combined[7], azimuth_radians, npts, transverse);
                vertical_basis(combined[0], combined[1], combined[2], azimuth_radians, npts, regional_basis[2]);
                /* east = R sin(az) + T cos(az), north = R cos(az) - T sin(az); overwrite R/T in place */
                for (k = 0; k < 5; k++) {
                    for (q = 0; q < npts; q++) {
                        float r = radial[k][q], t = transverse[k][q];
                        radial[k][q] = r * sin_az + t * cos_az;
                        transverse[k][q] = r * cos_az - t * sin_az;
                    }
                }
                for (c = 0; c < 3; c++) {
                    const int global_trace = p_count + sh_count + 3 * trace + c;
                    for (k = 0; k < 5; k++) {
                        float *source = regional_basis[c][k];
                        float *target = design_matrix + (size_t)(5 * subevent + k) * row_count + trace_start[global_trace];
                        for (q = 0; q < npts; q++) source[q] *= station->component_weight[c];
                        apply(source, npts, 1, setup.filter_sn[wave_regional], setup.filter_sd[wave_regional], setup.filter_sections[wave_regional]);
                        for (q = 0; q < npts; q++) target[q] = source[q] * setup.wave_weight[wave_regional];
                    }
                }
            }
        }
    }

    /* --- damped least squares with the total-moment-tensor rows (README section 4.3). The rows are
       weighted by (data norm / catalogue moment) so CMTscaling is dimensionless: a raw weight of 1e-3
       anchored an M8.8 but was invisible for an M6.9. */
    moment_cmt = sqrt((pow(setup.total_moment_tensor[0], 2) + pow(setup.total_moment_tensor[1], 2)
                       + pow(setup.total_moment_tensor[0] + setup.total_moment_tensor[1], 2)
                       + 2 * (pow(setup.total_moment_tensor[2], 2) + pow(setup.total_moment_tensor[3], 2) + pow(setup.total_moment_tensor[4], 2))) / 2);
    cmt_weight = setup.cmt_scaling * sqrt(data_energy) / (moment_cmt > 0 ? moment_cmt : 1);
    for (i = 0; i < unknown_count; i++) {
        const float *column_i = design_matrix + (size_t)i * row_count;
        for (k = 0; k <= i; k++) {
            const float *column_k = design_matrix + (size_t)k * row_count;
            double sum = 0;
            for (row = 0; row < row_count; row++) sum += (double)column_i[row] * column_k[row];
            if (i % 5 == k % 5) sum += cmt_weight * cmt_weight;
            if (i == k) sum += (double)setup.tikhonov_alpha * setup.tikhonov_alpha;
            normal_matrix[i * unknown_count + k] = normal_matrix[k * unknown_count + i] = sum;
        }
    }
    {
        /* relative ridge: co-located subevents with overlapping pulses make columns nearly collinear,
           and a subevent outside the windows leaves its columns dead; 1e-8 of the largest diagonal
           bounds every Cholesky pivot away from zero without visible bias */
        double largest_diagonal = 0;
        for (i = 0; i < unknown_count; i++) if (normal_matrix[i * unknown_count + i] > largest_diagonal) largest_diagonal = normal_matrix[i * unknown_count + i];
        for (i = 0; i < unknown_count; i++) normal_matrix[i * unknown_count + i] += 1e-8 * largest_diagonal;
    }
    if (!cholesky_factorize(normal_matrix, unknown_count)) goto finish;

    /* --- two linear solves: the second on data aligned with the lags of the first prediction */
    for (pass = 0; pass < 2; pass++) {
        const float *right_side_data = pass == 0 ? data_vector : shifted_data;
        for (i = 0; i < unknown_count; i++) {
            const float *column_i = design_matrix + (size_t)i * row_count;
            double sum = 0;
            for (row = 0; row < row_count; row++) sum += (double)column_i[row] * right_side_data[row];
            rhs[i] = sum + cmt_weight * cmt_weight * setup.total_moment_tensor[i % 5];
        }
        cholesky_solve(normal_matrix, unknown_count, rhs, solution);
        for (row = 0; row < row_count; row++) {
            double sum = 0;
            for (column = 0; column < unknown_count; column++) sum += solution[column] * design_matrix[(size_t)column * row_count + row];
            synthetic[row] = (float)sum;
        }
        if (pass == 1) break;
        for (trace = 0; trace < trace_total; trace++) {
            const int npts = trace_npts[trace], wave = trace_wave[trace];
            const float *d = data_vector + trace_start[trace], *s = synthetic + trace_start[trace];
            float *dx = shifted_data + trace_start[trace];
            int lag;
            float correlation;
            align_by_cross_correlation(d, s, npts, (int)(setup.max_shift_seconds[wave] / dt), &lag, &correlation);
            for (q = 0; q < npts; q++) dx[q] = (q - lag >= 0 && q - lag < npts) ? d[q - lag] : 0;
        }
    }
    for (i = 0; i < unknown_count; i++) moment_tensor[i] = solution[i];

    /* --- misfit (README section 4.4): residual energy after the final alignment, times
       exp(|1 - correlation|), normalized by the product of data and synthetic norms; the first and
       last 10 samples are left out (filter edges) */
    for (trace = 0; trace < trace_total; trace++) {
        const int npts = trace_npts[trace], wave = trace_wave[trace], edge = 10;
        const float *d = data_vector + trace_start[trace], *s = synthetic + trace_start[trace];
        double l2 = 0, data_norm = 0, synthetic_norm = 0;
        int lag;
        float correlation;
        align_by_cross_correlation(d, s, npts, (int)(setup.max_shift_seconds[wave] / dt), &lag, &correlation);
        for (q = edge; q < npts - edge; q++) {
            int shifted = q - lag;
            if (shifted < 0 || shifted >= npts) continue;
            l2 += ((double)d[shifted] - s[q]) * ((double)d[shifted] - s[q]);
            data_norm += (double)d[shifted] * d[shifted];
            synthetic_norm += (double)s[q] * s[q];
        }
        trace_lag[trace] = lag;
        trace_correlation[trace] = correlation;
        trace_l2[trace] = (float)l2;
        residual_sum[wave] += l2 * exp(fabs(1 - correlation));
        norm_sum[wave] += sqrt(data_norm * synthetic_norm);
    }

    /* --- penalties (README section 4.4): moment-weighted CLVD content, and the worst P/T-axis
       rotation of any subevent from the summed tensor beyond a 30 deg cone. Every subevent counts,
       however small, so a tiny subevent cannot absorb noise with an opposite mechanism; the 45 deg
       cone and moment-gated variants were tested in 2026-08 and rejected. */
    {
        double scalar_moment[max_subevents], clvd_fraction[max_subevents], strike[max_subevents][2], dip[max_subevents][2], rake[max_subevents][2];
        double total[6] = {0, 0, 0, 0, 0, 0}, total_clvd, total_strike[2], total_dip[2], total_rake[2];
        double moment_sum = 0, weighted_clvd = 0, worst_p = 0, worst_t = 0, p_total[3], t_total[3];
        double unit_moment = 1e27, data_misfit, stress_penalty;
        for (subevent = 0; subevent < subevent_count; subevent++) {
            double mxx = moment_tensor[5 * subevent], myy = moment_tensor[5 * subevent + 1], mxy = moment_tensor[5 * subevent + 2];
            double mxz = moment_tensor[5 * subevent + 3], myz = moment_tensor[5 * subevent + 4], mzz = -mxx - myy;
            scalar_moment[subevent] = sqrt((mxx * mxx + myy * myy + mzz * mzz + 2 * mxy * mxy + 2 * mxz * mxz + 2 * myz * myz) / 2);
            moment_sum += scalar_moment[subevent];
            mtdcmp_(&unit_moment, &mxx, &mxy, &mxz, &myy, &myz, &mzz, &clvd_fraction[subevent], strike[subevent], dip[subevent], rake[subevent]);
            total[0] += mxx; total[1] += mxy; total[2] += mxz; total[3] += myy; total[4] += myz; total[5] += mzz;
        }
        if (!(moment_sum > 0)) goto finish;
        mtdcmp_(&unit_moment, &total[0], &total[1], &total[2], &total[3], &total[4], &total[5], &total_clvd, total_strike, total_dip, total_rake);
        pressure_and_tension_axes(total_strike[0], total_dip[0], total_rake[0], p_total, t_total);
        for (subevent = 0; subevent < subevent_count; subevent++) {
            double p_axis[3], t_axis[3], p_angle, t_angle;
            weighted_clvd += sqrt(clvd_fraction[subevent]) * scalar_moment[subevent] / moment_sum;
            pressure_and_tension_axes(strike[subevent][0], dip[subevent][0], rake[subevent][0], p_axis, t_axis);
            p_angle = angle_between_axes_degrees(p_total, p_axis);
            t_angle = angle_between_axes_degrees(t_total, t_axis);
            if (p_angle > worst_p) worst_p = p_angle;
            if (t_angle > worst_t) worst_t = t_angle;
        }
        stress_penalty = (worst_p / 30 - 1 > 0 ? worst_p / 30 - 1 : 0) + (worst_t / 30 - 1 > 0 ? worst_t / 30 - 1 : 0);
        data_misfit = (residual_sum[0] + residual_sum[1] + residual_sum[2]) / (norm_sum[0] + norm_sum[1] + norm_sum[2]);
        *penalized_misfit = (float)(data_misfit * exp(weighted_clvd / setup.dc_constraint) * exp(stress_penalty));
        if (!is_finite_float(*penalized_misfit)) *penalized_misfit = 9999.0f;

        if (setup.outputs_enabled) {
            /* the weight-balancing step greps this line: shares of the residual sum per wave type */
            printf("bandresid: %e %e %e bandnorm: %e %e %e\n", residual_sum[0], residual_sum[1], residual_sum[2], norm_sum[0], norm_sum[1], norm_sum[2]);
            /* evaluate_misfit_components.sh greps this line: the three factors of the penalized misfit
               (data misfit, CLVD factor, stress factor), their product, and the worst P- and T-axis
               rotations in degrees (unlabeled until 2026-11-16) */
            printf("misfit_components: %f %f %f %f %f %f\n", data_misfit, exp(weighted_clvd / setup.dc_constraint), exp(stress_penalty), *penalized_misfit, worst_p, worst_t);
            fflush(stdout);
            if (!getenv("FFWD_COMPACT")) {
                const char *wave_header[3] = {"# P waves", "# SH waves", "# Rayleigh waves"};
                const char *obs_prefix[5] = {"P_obs_", "SH_obs_", "rayle_obs_", "rayln_obs_", "raylz_obs_"};
                const char *syn_prefix[5] = {"P_syn_", "SH_syn_", "rayle_syn_", "rayln_syn_", "raylz_syn_"};
                FILE *moment_file = fopen("fm.dat", "w"), *residual_file = fopen("sta_residual.dat", "w");
                float *aligned = malloc(sizeof(float) * (max_npts + 1));
                float *pulse = calloc(setup.stf_npts + 1, sizeof(float)), *pulse_sum = calloc(setup.stf_npts + 1, sizeof(float));
                char file_name[128];
                int wave, current_wave = -1;
                if (moment_file == NULL || residual_file == NULL) fail("cannot write fm.dat / sta_residual.dat", NULL);
                remove_sac_files("stf");
                remove_sac_files("waveforms");
                for (subevent = 0; subevent < subevent_count; subevent++) {
                    const double *m = moment_tensor + 5 * subevent;
                    const float duration = model[7 * subevent + 3], centre = model[7 * subevent];
                    const float sigma = duration / (2 * (float)sqrt(2 * log(10)));
                    printf("1e27 %.4f %.4f %.4f %.4f %.4f %.4f \n", m[0], m[2], m[3], m[1], m[4], -m[0] - m[1]);
                    fprintf(moment_file, "1e27 %.4f %.4f %.4f %.4f %.4f %.4f\n", m[0], m[2], m[3], m[1], m[4], -m[0] - m[1]);
                    /* stf/stf_NNNN.sac: unit-area pulse per subevent; mrf_sum.sac: moment-rate function */
                    for (q = 0; q < setup.stf_npts; q++) {
                        float time = setup.stf_begin + q * dt;
                        float value = 1 / ((float)sqrt(2 * 3.1416) * sigma) * expf(-(time - centre) * (time - centre) / (2 * sigma * sigma));
                        pulse[q] = value * dt;
                        pulse_sum[q] += value * (float)scalar_moment[subevent];
                    }
                    snprintf(file_name, sizeof file_name, "stf/stf_%04d.sac", subevent);
                    write_sac_trace(file_name, pulse, setup.stf_npts);
                    printf("clvd[%d] = %f\n", subevent + 1, clvd_fraction[subevent]);
                }
                printf("\n1e27 %.4f %.4f %.4f %.4f %.4f %.4f \n", total[0], total[1], total[2], total[3], total[4], total[5]);
                write_sac_trace("stf/mrf_sum.sac", pulse_sum, setup.stf_npts);
                for (trace = 0; trace < trace_total; trace++) {
                    const int npts = trace_npts[trace], lag = trace_lag[trace];
                    const float *d = data_vector + trace_start[trace], *s = synthetic + trace_start[trace];
                    int local, kind;
                    wave = trace_wave[trace];
                    if (wave != current_wave) { fprintf(residual_file, "%s\n", wave_header[wave]); current_wave = wave; }
                    if (wave == wave_p) { local = trace; kind = 0; }
                    else if (wave == wave_sh) { local = trace - p_count; kind = 1; }
                    else { local = trace - p_count - sh_count; kind = 2 + local % 3; local /= 3; }
                    for (q = 0; q < npts; q++) aligned[q] = (q - lag >= 0 && q - lag < npts) ? d[q - lag] : 0;
                    snprintf(file_name, sizeof file_name, "waveforms/%s%04d.sac", obs_prefix[kind], local);
                    write_sac_trace(file_name, aligned, npts);
                    snprintf(file_name, sizeof file_name, "waveforms/%s%04d.sac", syn_prefix[kind], local);
                    write_sac_trace(file_name, s, npts);
                    {
                        /* per-trace share of the residual, scaled so that 1 = average trace of its type */
                        int type_count = wave == wave_p ? p_count : (wave == wave_sh ? sh_count : 3 * regional_count);
                        int index_in_type = wave == wave_p ? trace : (wave == wave_sh ? trace - p_count : trace - p_count - sh_count);
                        double share = norm_sum[wave] > 0 ? trace_l2[trace] / norm_sum[wave] * exp(fabs(1 - trace_correlation[trace])) * type_count : 0;
                        fprintf(residual_file, "%d %f %f %f\n", index_in_type + 1, share, trace_correlation[trace], lag * dt);
                    }
                }
                fclose(moment_file);
                fclose(residual_file);
                free(aligned); free(pulse); free(pulse_sum);
                {
                    double norm_total = norm_sum[0] + norm_sum[1] + norm_sum[2];
                    printf("residualp= %f  residualsh= %f residualrayl= %f\n", residual_sum[0] / norm_total, residual_sum[1] / norm_total, residual_sum[2] / norm_total);
                    printf("dataresidual= %f  finalresidual= %f\n", data_misfit, *penalized_misfit);
                }
            }
        }
    }
finish:
    free(trace_npts); free(trace_start); free(trace_wave); free(data_vector); free(shifted_data); free(synthetic);
    free(design_matrix); free(moment_rate); free(work); free(trace_correlation); free(trace_l2); free(trace_lag);
}
