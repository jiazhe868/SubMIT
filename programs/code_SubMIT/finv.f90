! finv.f90 - Bayesian MCMC sampling of N subevents (time, position, depth, duration) with MPI chains.
!
! What it does. Every MPI rank runs one Markov chain over the nonlinear parameters of a fixed number
! of subevents (the moment tensors are solved exactly inside forward_model.c). The posterior is
! p(theta|d) ~ exp(-phi(theta) / (2 sigma^2 phi_ref)) * prior(theta) with sigma = 0.05 (README 4.5);
! the prior is the aftershock density map (seisdens.dat) times hard rules: search bounds, centroid
! times increasing by at least t_gap, no moment before the origin (t_i >= T_i/2), causality
! (|x_i - x_1| <= vmax_cau (t_i - s_1), with s_1 = t_1 - T_1/2 the start of E1) and the 45-degree
! corner rule for three or more subevents.
! Each chain writes its post-burn-in samples (NNNchain.dat) and its 50 best samples (NNNbest.dat).
!
! The first subevent is where the rupture starts (two or more subevents): it sits at the hypocentre,
! epicentre and depth, and the rupture front that the causality rule limits leaves it when its moment
! rate begins. Until 2026-11-17 E1 was pinned at the epicentre only, its depth was searched, and
! causality was measured from E1's centroid, |x_i - x_1| <= vmax_cau (t_i - t_1), with a cap of
! 1.5 Vs for shallow strike-slip events. Measured on the v3 runs of the examples (PITFALLS D6):
! - Chile (hypocentre at 122.9 km): E1 at 137-146 km in every run of two seeds. With the depth pinned
!   and causality still measured from E1's centroid, the 4-subevent misfit rose by 26% (0.9207
!   against 0.7318) and the fourth subevent went to 75 km at 23 s: the deepest moment (168 km at 10 s)
!   cannot be reached from the centroid of E1 (2-4 s) at the S speed, but it can from its start.
! - Venezuela: E1 at the epicentre at 15.8 s, E2 127 km away at 38.7 s, 5.5 km/s, inside the 1.5 Vs
!   (5.85 km/s) that the cap allowed. That E1 was a 4.5-s pulse of Mw 6.5 worth 0.1% of misfit: with
!   E1 early the misfit was 2.0591 against 2.0562.
! A stricter rule was tried and dropped: E1 starting at the origin time (t_1 = T_1/2), causality from
! the origin time. It fits Chile and Venezuela as well, but an earthquake whose first resolvable
! moment release starts late needs one more subevent for the silent beginning: Colombia (17 s) fit
! 1.6229 with two subevents and 1.6411 with three, against 1.4660 and 1.3989 with the rule above, and
! Kamchatka (18 s) kept six subevents for the part that seven had described (0.9065 against 0.8333
! at the limit of seven).
!
! Why a population of chains. The penalized misfit is extremely peaked (a model 1% worse has e^-2
! of the density) and multimodal, and independent chains settle in different local modes: in the
! published runs typically only 2-6 of 32 chains came within 2% of the best misfit, the median
! chain's best was 4-13% worse, and the split R-hat of the gated chains reached 2-11 (Chile 4sub
! E1 depth 2.1, California 3sub E3 x 11.0). The best misfit per N, and with it the selected number
! of subevents, then depended on luck: the 2026 Colombia runs selected 5, 4, 3, 2, 4 subevents under
! small configuration changes, and in the 2026 Panama test N = 2 rested on one chain of 16.
! The sampler now has five parts:
! - the burn-in is split into anneal_stages stages; between stages every chain copies the state of a
!   chain drawn with weight exp(-dphi / (2 sigma^2 phi_ref)), its posterior density, among the chains
!   of its island of 8; the last round draws from all chains. Chains in poor modes are replaced by
!   copies of good ones, and islands keep several modes alive until that last round;
! - 5% of the proposals move one subevent to a position drawn from the prior (relocation), so a
!   chain can jump between modes tens of km apart that local moves cannot connect;
! - half of the chains start from the best N-1 model plus one subevent (warm start,
!   warm_start.dat), so the search for N begins from the model already found with N-1. The new
!   subevent is drawn from the prior and needs a free slot in time (2.5 s from its neighbours).
!   When no two subevents of the N-1 model leave such a slot between them, half of the draws
!   insert it between two neighbours and move the centroid times of the others apart
!   (draw_warm_start_model, with the Chile measurement). Before the insertion existed, a model
!   without a slot started every chain of the next N from the prior: the 4-subevent model of an
!   early 32-rank California run (centroid times 5.2, 8.5, 12.8 and 15.4 s in a range that ended
!   at 16 s) did so for its 5- and 6-subevent runs. The log line of each chain says which start
!   it got;
! - island_best_misfits.txt records each island's best before the last round: the islands are
!   independent searches until then, and plot_lcurve.py lists their bests next to the L-curve;
! - a random-walk proposal that breaks a hard rule is a rejection (PITFALLS F1).
! Measured on the California 3-subevent configuration (16 ranks, 9000 + 3000 iterations, seeds 0
! and 500): independent chains (3000-iteration burn-in) left 1 of 16 chains within 2% of the best
! and the median chain 34% worse; with resampling 10-16 of 16 chains ended within 2% in seven of
! eight runs. The final design reached 1.840 and 1.844 (seeds 0 and 500, 878-916 s); the version
! that redrew infeasible proposals reached 1.761, 1.840 and 1.847, and a single population cooled
! from anneal_t0 = 3 reached 1.871 and 1.933. A mode 4.5% below the usual one (E3 east of the
! epicentre instead of west) was found in 1 of 5 runs, close to the 5% selection band, which is
! why the island bests are recorded instead of trusting the agreement of the chains.
! Sampling after the burn-in is plain Metropolis-Hastings at temperature 1.
!
! History of measured defects fixed here (v2 -> v3):
! - the single-parameter moves picked the subevent with 2+ANINT(u*(neq-2)): for one subevent half of
!   the x/y moves changed an inactive slot (wasted forward runs), for three or more the first and
!   last subevents were proposed half as often as the others. Now uniform.
! - the prior draw indexed the nonzero density cells from 0 (one cell out of bounds).
! - for one subevent the adaptive (joint) moves and the burn-in consolidation ignored x and y of E1,
!   and the consolidation copied a peer's misfit onto a state with this chain's own x and y: the
!   recorded misfit did not belong to the recorded model (Calama 1sub y: R-hat 4.7).
! - the one-shot consolidation through MPI one-sided windows (10% band) hung at MPI_Win_free on the
!   Kamchatka run; collectives at fixed iteration counts replace it.
! - a NaN/Inf misfit or a zero prior density could pass the acceptance test under -Ofast.
! - output units equal to the rank number collided with stdin/stdout for ranks 5 and 6.
!
! Run (step3 writes runmpi.sh):  make finv && mpirun -np 32 ./finv > mpi.out
! search_par.file keys: burn_in nsample neq_min neq_max x_min x_max y_min y_max z_min z_max
!   dura_min dura_max vr_min vr_max theta_min theta_max cen_min cen_max vmax_cau
!   optional: anneal_stages (default 4: population resampled at 25, 50, 75% of the burn-in;
!   1 = independent chains), anneal_t0 (default 1: every stage at temperature 1; > 1 = stages
!   cooled geometrically from anneal_t0 to 1)
program finv
    implicit none
    include 'mpif.h'
    integer, parameter :: max_subevents = 9
    integer, parameter :: max_adaptive = 5 * max_subevents      ! active parameters of the adaptive moves
    integer, parameter :: best_samples_written = 50
    integer, parameter :: display_interval = 250
    integer, parameter :: adaptive_start = 500                   ! states before the adaptive moves begin
    integer, parameter :: island_size = 8                        ! chains per island of the population resampling
    ! posterior width (README 4.5); plot_lcurve.py uses the same value in the Akaike criterion that
    ! selects the number of subevents, and tests/test_regressions.py checks that the two are equal
    real, parameter :: sigma_relative = 0.05
    real, parameter :: t_gap_nominal = 2.5                       ! minimum gap of centroid times (s)
    real, parameter :: degrees_per_radian = 57.29578

    ! --- configuration (search_par.file)
    integer :: burn_in, nsample, neq_min, neq_max, anneal_stages
    real :: x_min, x_max, y_min, y_max, z_min, z_max, dura_min, dura_max, vr_min, vr_max
    real :: theta_min, theta_max, cen_min, cen_max, vmax_cau, anneal_t0
    ! --- aftershock prior on the regular grid of seisdens.dat
    integer :: grid_cell_count, nonzero_cell_count
    real, allocatable :: grid_x(:), grid_y(:), grid_density(:), nonzero_x(:), nonzero_y(:)
    real :: grid_x_first, grid_x_last, grid_y_first, grid_y_last, grid_spacing
    ! --- chain state: para(i, :) = x, y, z, duration, rupture velocity, rupture direction, centroid time
    real :: para(max_subevents, 7), para_proposed(max_subevents, 7), para_best(max_subevents, 7)
    real :: like, like_proposed, like_best, like_reference, temperature
    real :: model_vector(7 * max_subevents), source_s_velocity
    ! --- hypocentre: the depth that the arrival picks of the data refer to (evdp of the headers), moved
    ! into the depth range; a negative header depth (hand-made data without evdp) leaves E1's depth free
    real :: hypocentre_header_depth, hypocentre_depth
    logical :: first_depth_is_pinned
    integer :: neq
    ! --- proposal widths and adaptive-Metropolis state (Haario et al. 2001)
    real :: width_x, width_y, width_z, width_duration, width_centroid
    integer :: adaptive_count, adaptive_subevent(max_adaptive), adaptive_column(max_adaptive), adaptive_ready
    real :: adaptive_width(max_adaptive), adaptive_mean(max_adaptive), adaptive_m2(max_adaptive, max_adaptive)
    real :: adaptive_factor(max_adaptive, max_adaptive), adaptive_state_count
    real :: joint_log_scale, direction_log_scale, adaptive_vector(max_adaptive), adaptive_delta(max_adaptive)
    real :: gaussian_draws(max_adaptive), direction(max_adaptive), step_length
    integer :: joint_proposed, joint_accepted, direction_proposed, direction_accepted
    integer :: relocation_proposed, relocation_accepted
    integer :: single_proposed(7), single_accepted(7)
    ! --- population annealing
    integer :: stage, stage_length, state_size
    real, allocatable :: population(:, :), own_state(:)
    ! --- bookkeeping
    integer :: rank, rank_count, mpi_status_code, seed_offset, random_state
    integer :: iteration, sample_count, redraw_count, infeasible_count, move_kind, subevent_index
    integer :: i, j, k, read_status, chain_unit, best_unit
    real :: u, t_gap, centroid_bin_width, prior_current, prior_proposed, log_acceptance, time_started, time_finished
    real :: initial_draw_like, temporary(3)
    ! --- warm start: the best model of the run with one subevent fewer (warm_start.dat, written by step 3)
    integer :: warm_subevent_count
    real :: warm_para(max_subevents, 7)
    logical :: is_warm_chain, warm_draw_ok, started_warm
    real, allocatable :: samples(:, :, :)
    real :: sample_swap(max_subevents, 8)
    logical :: is_feasible, accepted
    character(len=20) :: parameter_name
    character(len=16) :: rank_text
    real :: parameter_value
    real, external :: gaussian_deviate, uniform_deviate, seismicity_prior

    ! ---------------------------------------------------------------------------------------------
    ! configuration and prior
    ! ---------------------------------------------------------------------------------------------
    vmax_cau = -1
    anneal_t0 = 1
    anneal_stages = 4
    open(10, file='search_par.file', status='old')
    do
        read(10, *, iostat=read_status) parameter_name, parameter_value
        if (read_status /= 0) exit
        select case (trim(parameter_name))
        case ('burn_in');   burn_in = int(parameter_value)
        case ('nsample');   nsample = int(parameter_value)
        case ('neq_min');   neq_min = int(parameter_value)
        case ('neq_max');   neq_max = int(parameter_value)
        case ('x_min');     x_min = parameter_value
        case ('x_max');     x_max = parameter_value
        case ('y_min');     y_min = parameter_value
        case ('y_max');     y_max = parameter_value
        case ('z_min');     z_min = parameter_value
        case ('z_max');     z_max = parameter_value
        case ('dura_min');  dura_min = parameter_value
        case ('dura_max');  dura_max = parameter_value
        case ('vr_min');    vr_min = parameter_value
        case ('vr_max');    vr_max = parameter_value
        case ('theta_min'); theta_min = parameter_value
        case ('theta_max'); theta_max = parameter_value
        case ('cen_min');   cen_min = parameter_value
        case ('cen_max');   cen_max = parameter_value
        case ('vmax_cau');  vmax_cau = parameter_value
        case ('anneal_t0'); anneal_t0 = max(1.0, parameter_value)
        case ('anneal_stages'); anneal_stages = max(1, int(parameter_value))
        case default;       print *, 'Warning: unknown search_par.file key ', trim(parameter_name)
        end select
    end do
    close(10)
    if (neq_min /= neq_max) stop 'finv: neq_min must equal neq_max (one subevent count per run)'
    if (neq_min < 1 .or. neq_min > max_subevents) stop 'finv: neq must be 1..9'
    neq = neq_min
    open(13, file='rank.dat', status='old')
    read(13, *) seed_offset
    close(13)
    call read_warm_start()

    ! seisdens.dat: x y density on a regular grid (x fastest), written by caldens.py
    grid_cell_count = 0
    nonzero_cell_count = 0
    open(12, file='seisdens.dat', status='old')
    do
        read(12, *, iostat=read_status) temporary
        if (read_status /= 0) exit
        grid_cell_count = grid_cell_count + 1
        if (temporary(3) > 0) nonzero_cell_count = nonzero_cell_count + 1
    end do
    if (nonzero_cell_count == 0) stop 'finv: seisdens.dat has no cell with positive density'
    allocate(grid_x(grid_cell_count), grid_y(grid_cell_count), grid_density(grid_cell_count))
    allocate(nonzero_x(nonzero_cell_count), nonzero_y(nonzero_cell_count))
    rewind(12)
    j = 0
    do i = 1, grid_cell_count
        read(12, *) grid_x(i), grid_y(i), grid_density(i)
        if (grid_density(i) > 0) then
            j = j + 1
            nonzero_x(j) = grid_x(i)
            nonzero_y(j) = grid_y(i)
        end if
    end do
    close(12)
    open(11, file='seisgrids.dat', status='old')
    read(11, *) grid_x_first, grid_x_last, grid_y_first, grid_y_last, grid_spacing
    close(11)

    ! ---------------------------------------------------------------------------------------------
    ! MPI, data, proposal set-up
    ! ---------------------------------------------------------------------------------------------
    call cpu_time(time_started)
    call MPI_INIT(mpi_status_code)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, rank_count, mpi_status_code)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_status_code)
    ! one random stream per chain, seeded by rank + rank.dat offset (a different offset = a rerun
    ! with independent chains)
    random_state = rank + seed_offset
    call sub_init(source_s_velocity, hypocentre_header_depth)
    first_depth_is_pinned = hypocentre_header_depth >= 0
    hypocentre_depth = max(z_min, min(z_max, hypocentre_header_depth))
    if (rank == 0 .and. first_depth_is_pinned) write(*, '(A, F8.2, A, F8.2, A)') ' hypocentre depth of the data headers ', &
        hypocentre_header_depth, ' km; E1 is placed at ', hypocentre_depth, ' km'
    if (rank == 0 .and. .not. first_depth_is_pinned) write(*, *) 'WARNING: no hypocentre depth in the data headers: ', &
        'the depth of E1 is searched'
    if (vmax_cau < 0) vmax_cau = abs(vmax_cau) * source_s_velocity     ! negative = factor of Vs
    if (rank == 0) write(*, *) 'causality speed cap vmax_cau (km/s):', vmax_cau

    ! short centroid ranges: shrink the gap so the staggered initial bins keep a positive width
    t_gap = t_gap_nominal
    centroid_bin_width = ((cen_max - cen_min) - (neq - 1) * t_gap) / neq
    if (centroid_bin_width <= 0) then
        t_gap = 0.5 * (cen_max - cen_min) / max(1, neq - 1)
        centroid_bin_width = ((cen_max - cen_min) - (neq - 1) * t_gap) / neq
        if (rank == 0) write(*, *) 'WARNING: centroid range too short for', neq, 'subevents; gap reduced to', t_gap
    end if

    ! proposal widths scale with the search box (fixed widths mixed far too slowly on +/-100 km boxes);
    ! the coefficients reproduce the historical 2/3/1/0.5/0.5 at Calama-scale ranges
    width_x = max(0.033 * (x_max - x_min), 0.5)
    width_y = max(0.040 * (y_max - y_min), 0.5)
    width_z = max(0.025 * (z_max - z_min), 0.5)
    width_duration = max(0.030 * (dura_max - dura_min), 0.25)
    width_centroid = max(0.020 * (cen_max - cen_min), 0.25)
    if (rank == 0) write(*, *) 'proposal widths x,y,z,dur,cen:', width_x, width_y, width_z, width_duration, width_centroid

    ! active parameters of the adaptive moves: x, y, depth, duration, centroid time; not the position
    ! of E1 when it is pinned at the hypocentre, i.e. for two or more subevents
    adaptive_count = 0
    do i = 1, neq
        do k = 1, 5
            j = k
            if (k == 4) j = 4
            if (k == 5) j = 7
            if (neq > 1 .and. i == 1 .and. k <= 2) cycle
            if (neq > 1 .and. i == 1 .and. k == 3 .and. first_depth_is_pinned) cycle
            adaptive_count = adaptive_count + 1
            adaptive_subevent(adaptive_count) = i
            adaptive_column(adaptive_count) = j
            select case (j)
            case (1); adaptive_width(adaptive_count) = width_x
            case (2); adaptive_width(adaptive_count) = width_y
            case (3); adaptive_width(adaptive_count) = width_z
            case (4); adaptive_width(adaptive_count) = width_duration
            case (7); adaptive_width(adaptive_count) = width_centroid
            end select
        end do
    end do
    adaptive_state_count = 0
    adaptive_mean = 0
    adaptive_m2 = 0
    adaptive_ready = 0
    joint_log_scale = 0
    direction_log_scale = 0
    joint_proposed = 0; joint_accepted = 0; direction_proposed = 0; direction_accepted = 0
    relocation_proposed = 0; relocation_accepted = 0
    single_proposed = 0; single_accepted = 0
    ! population state: misfit, the 7 parameters of every subevent, the adaptive statistics
    state_size = 1 + 7 * max_subevents + 4 + 2 * max_adaptive + 2 * max_adaptive * max_adaptive
    allocate(population(state_size, rank_count), own_state(state_size))

    ! ---------------------------------------------------------------------------------------------
    ! initial model: best of 101 feasible draws from the prior
    ! ---------------------------------------------------------------------------------------------
    like = 1e30
    ! (warm_subevent_count >= 1: without warm_start.dat the count is 0, which equals neq - 1 for one
    ! subevent and once made the odd chains of every N = 1 run start at the epicentre)
    is_warm_chain = warm_subevent_count >= 1 .and. warm_subevent_count == neq - 1 .and. mod(rank, 2) == 1
    started_warm = .false.
    do j = 1, 101
        warm_draw_ok = .false.
        if (is_warm_chain) call draw_warm_start_model(para_proposed, warm_draw_ok)
        if (.not. warm_draw_ok) call draw_initial_model(para_proposed)
        call model_to_vector(neq, para_proposed, model_vector)
        call sub_forward(neq, model_vector, initial_draw_like)
        if (initial_draw_like < like) then
            like = initial_draw_like
            para = para_proposed
            started_warm = warm_draw_ok
        end if
    end do
    if (like >= 9999.0) write(*, *) 'WARNING rank', rank + seed_offset, ': no initial draw had a finite misfit'
    like_best = like
    para_best = para
    ! "warm" only when the start is a warm draw: a warm chain whose N-1 model leaves no room for
    ! another subevent starts from the prior like the others, and was once labelled warm all the same
    write(*, *) 'Rank ', rank + seed_offset, 'Like_init', like, merge(' warm', ' cold', started_warm)
    call MPI_ALLREDUCE(like, like_reference, 1, MPI_REAL, MPI_MIN, MPI_COMM_WORLD, mpi_status_code)

    ! ---------------------------------------------------------------------------------------------
    ! sampling: annealed burn-in with population resampling, then Metropolis-Hastings at T = 1
    ! ---------------------------------------------------------------------------------------------
    if (anneal_stages == 1) like_reference = like_best
    stage_length = max(1, burn_in / anneal_stages)
    allocate(samples(nsample, neq, 8))
    samples = 0
    sample_count = 0
    infeasible_count = 0
    stage = 0
    temperature = stage_temperature(0)
    do iteration = 1, burn_in + nsample
        ! stage boundary: resample the population, then continue at the next (lower) temperature
        if (iteration <= burn_in .and. stage < anneal_stages - 1 .and. iteration == (stage + 1) * stage_length + 1) then
            ! islands of 8 chains resample among themselves, the last round across all chains
            call resample_population(stage_temperature(stage), stage_temperature(stage + 1), stage + 1 == anneal_stages - 1)
            stage = stage + 1
            temperature = stage_temperature(stage)
            if (rank == 0) write(*, *) 'annealing stage', stage + 1, 'of', anneal_stages, 'temperature', temperature, &
                                       'reference misfit', like_reference
        end if
        if (iteration == burn_in + 1) temperature = 1.0

        ! --- propose. A random-walk move (kinds 1-9) is proposed once, and a proposal that breaks a
        ! hard rule is a rejection, as Metropolis-Hastings requires. v2 redrew until every rule held,
        ! which samples pi(x) Z(x) instead of pi(x) (Z = the share of proposals from x that are
        ! feasible) and under-weights states near the bounds and the time-gap and causality rules. A
        ! relocation (kind 10) draws the subevent's position from the prior, independently of where it
        ! is; redrawing it until the rules hold samples the prior restricted to a set fixed by the other
        ! subevents, which the move leaves unchanged, so that redraw stays exact.
        call propose(move_kind, subevent_index)
        call place_first_subevent(para_proposed)
        is_feasible = satisfies_hard_rules(para_proposed)
        redraw_count = 0
        do while (move_kind == 10 .and. .not. is_feasible .and. redraw_count < 200)
            call draw_relocation(subevent_index)
            is_feasible = satisfies_hard_rules(para_proposed)
            redraw_count = redraw_count + 1
        end do

        ! --- evaluate and accept (prior ratio times likelihood ratio at the stage temperature); an
        ! infeasible proposal needs no forward run
        accepted = .false.
        if (is_feasible) then
            call model_to_vector(neq, para_proposed, model_vector)
            call sub_forward(neq, model_vector, like_proposed)
            if (.not. is_finite_real(like_proposed)) like_proposed = 9999.0
            prior_current = 1
            prior_proposed = 1
            do i = 1, neq
                ! E1 is pinned at the hypocentre for two or more subevents: its density is a constant
                ! factor that cancels in the ratio, and where the aftershock map is zero at the epicentre
                ! it would make every state impossible and freeze all chains at their initial draw
                if (neq > 1 .and. i == 1) cycle
                prior_current = prior_current * seismicity_prior(para(i, 1), para(i, 2), grid_spacing, grid_cell_count, &
                                                                 grid_x, grid_y, grid_density)
                prior_proposed = prior_proposed * seismicity_prior(para_proposed(i, 1), para_proposed(i, 2), grid_spacing, &
                                                                   grid_cell_count, grid_x, grid_y, grid_density)
            end do
            if (prior_proposed <= 0) then
                accepted = .false.                     ! zero prior density: never
            else if (prior_current <= 0) then
                accepted = .true.                      ! leaving a zero-density start
            else
                log_acceptance = (like - like_proposed) / (2 * sigma_relative**2 * like_reference * temperature)
                if (move_kind /= 10) log_acceptance = log_acceptance + log(prior_proposed) - log(prior_current)
                accepted = log(max(uniform_deviate(random_state), 1e-30)) < log_acceptance
            end if
        else
            infeasible_count = infeasible_count + 1
        end if
        if (move_kind == 10) then
            relocation_proposed = relocation_proposed + 1
            if (accepted) relocation_accepted = relocation_accepted + 1
        end if
        if (move_kind >= 1 .and. move_kind <= 7) then
            single_proposed(move_kind) = single_proposed(move_kind) + 1
            if (accepted) single_accepted(move_kind) = single_accepted(move_kind) + 1
        end if
        if (move_kind == 8) then
            joint_proposed = joint_proposed + 1
            if (accepted) joint_accepted = joint_accepted + 1
            call adapt_scale(joint_log_scale, joint_proposed, merge(1.0, 0.0, accepted), 0.234)
        end if
        if (move_kind == 9) then
            direction_proposed = direction_proposed + 1
            if (accepted) direction_accepted = direction_accepted + 1
            call adapt_scale(direction_log_scale, direction_proposed, merge(1.0, 0.0, accepted), 0.44)
        end if
        if (accepted) then
            para = para_proposed
            like = like_proposed
            if (like < like_best) then
                like_best = like
                para_best = para
                if (anneal_stages == 1) like_reference = like_best
            end if
        end if
        call update_adaptive_statistics()

        ! --- store post-burn-in samples
        if (iteration > burn_in) then
            sample_count = sample_count + 1
            do i = 1, neq
                samples(sample_count, i, :) = (/ like, para(i, 7), para(i, 1), para(i, 2), para(i, 4), para(i, 5), &
                                                 para(i, 6), para(i, 3) /)
            end do
        end if
        if (mod(iteration, display_interval) == 0) call print_progress()
    end do

    ! ---------------------------------------------------------------------------------------------
    ! outputs: all post-burn-in samples, and the 50 best
    ! ---------------------------------------------------------------------------------------------
    write(rank_text, '(I0)') rank + seed_offset
    chain_unit = 21
    best_unit = 22
    open(chain_unit, file=trim(rank_text)//'chain.dat', status='replace')
    do j = 1, sample_count
        call write_sample_line(chain_unit, j)
    end do
    close(chain_unit)
    ! partial selection sort: only the best 50 are needed
    do i = 1, min(best_samples_written, sample_count)
        k = i
        do j = i + 1, sample_count
            if (samples(j, 1, 1) < samples(k, 1, 1)) k = j
        end do
        if (k /= i) then
            sample_swap(1:neq, :) = samples(i, :, :)
            samples(i, :, :) = samples(k, :, :)
            samples(k, :, :) = sample_swap(1:neq, :)
        end if
    end do
    open(best_unit, file=trim(rank_text)//'best.dat', status='replace')
    do j = 1, min(best_samples_written, sample_count)
        call write_sample_line(best_unit, j)
    end do
    close(best_unit)
    call cpu_time(time_finished)
    if (rank == 0) write(*, *) 'time taken by the code was', time_finished - time_started, 'seconds'
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_status_code)
    call MPI_FINALIZE(mpi_status_code)

contains

    ! temperature of annealing stage s (0-based): anneal_t0 at the start, 1 in the last stage
    real function stage_temperature(stage_index)
        integer, intent(in) :: stage_index
        if (anneal_stages <= 1) then
            stage_temperature = 1.0
        else
            stage_temperature = anneal_t0 ** (real(anneal_stages - 1 - stage_index) / real(anneal_stages - 1))
        end if
    end function

    ! one draw from the prior that satisfies every hard rule: positions from the density map,
    ! depths uniform, centroid times in staggered bins, E1 at the hypocentre for two or more subevents
    subroutine draw_initial_model(new_para)
        real, intent(out) :: new_para(max_subevents, 7)
        integer :: attempt, cell, n, nearest, position_draw
        real :: x_draw, y_draw, duration_deviate, distance(max_subevents), swap(3), separation, rupture_start
        do attempt = 1, 100000
            new_para = 0
            do n = 1, neq
                ! rejection sampling of a position from the density map; cells are indexed 1..count
                do position_draw = 1, 100
                    cell = 1 + min(nonzero_cell_count - 1, int(uniform_deviate(random_state) * nonzero_cell_count))
                    x_draw = nonzero_x(cell) + (uniform_deviate(random_state) - 0.5) * grid_spacing
                    y_draw = nonzero_y(cell) + (uniform_deviate(random_state) - 0.5) * grid_spacing
                    if (uniform_deviate(random_state) < seismicity_prior(x_draw, y_draw, grid_spacing, grid_cell_count, &
                                                                       grid_x, grid_y, grid_density)) exit
                end do
                new_para(n, 1) = max(min(x_draw, x_max), x_min)
                new_para(n, 2) = max(min(y_draw, y_max), y_min)
                new_para(n, 3) = z_min + uniform_deviate(random_state) * (z_max - z_min)
                duration_deviate = uniform_deviate(random_state)
                new_para(n, 5) = 0.01                     ! point sources (README 4.1)
                new_para(n, 6) = 0
                new_para(n, 7) = cen_min + (n - 1) * (centroid_bin_width + t_gap) &
                                 + uniform_deviate(random_state) * centroid_bin_width
                new_para(n, 4) = dura_min + duration_deviate * (min(new_para(n, 7) * 2, dura_max) - dura_min)
            end do
            if (neq > 1) then
                call place_first_subevent(new_para)      ! E1 pinned at the hypocentre
                ! label subevents by distance from E1 (nearest first), positions only
                do n = 2, neq
                    distance(n) = sqrt(sum((new_para(n, 1:3) - new_para(1, 1:3))**2))
                end do
                do n = 2, neq - 1
                    nearest = n - 1 + minloc(distance(n:neq), 1)
                    if (nearest /= n) then
                        swap = new_para(n, 1:3); new_para(n, 1:3) = new_para(nearest, 1:3); new_para(nearest, 1:3) = swap
                        separation = distance(n); distance(n) = distance(nearest); distance(nearest) = separation
                    end if
                end do
                ! a start farther from E1 than causality allows is unreachable by single moves: delay it,
                ! or pull it toward E1
                rupture_start = rupture_start_time(new_para)
                do n = 2, neq
                    separation = sqrt(sum((new_para(n, 1:3) - new_para(1, 1:3))**2))
                    if (separation > vmax_cau * (new_para(n, 7) - rupture_start)) then
                        new_para(n, 7) = min(cen_max, rupture_start + separation / vmax_cau + 0.2)
                        if (separation > vmax_cau * (new_para(n, 7) - rupture_start)) then
                            new_para(n, 1:3) = new_para(1, 1:3) + (new_para(n, 1:3) - new_para(1, 1:3)) &
                                               * 0.9 * vmax_cau * (new_para(n, 7) - rupture_start) / separation
                            new_para(n, 3) = max(min(new_para(n, 3), z_max), z_min)
                        end if
                    end if
                end do
            end if
            if (satisfies_hard_rules(new_para)) return
        end do
        stop 'finv: no feasible initial model in 100000 draws (search bounds too tight?)'
    end subroutine

    ! warm_start.dat: one line of finv's sample format (best_model_hybrid.dat of the N-1 run)
    subroutine read_warm_start()
        integer :: status, n, subevent
        real :: values(1 + 8 * max_subevents)
        character(len=4096) :: line
        warm_subevent_count = 0
        open(14, file='warm_start.dat', status='old', iostat=status)
        if (status /= 0) return
        read(14, '(A)', iostat=status) line
        close(14)
        if (status /= 0 .or. index(line, '<neq') == 0) return
        read(line(1:index(line, '<neq') - 1), *, iostat=status) n
        if (status /= 0 .or. n < 1 .or. n >= max_subevents) return
        read(line(index(line, '<neq') + 4:), *, iostat=status) values(1:8 * n)
        if (status /= 0) return
        warm_para = 0
        do subevent = 1, n
            ! file order per subevent: misfit cen x y dur vr theta z
            warm_para(subevent, 1) = values(8 * (subevent - 1) + 3)
            warm_para(subevent, 2) = values(8 * (subevent - 1) + 4)
            warm_para(subevent, 3) = values(8 * (subevent - 1) + 8)
            warm_para(subevent, 4) = values(8 * (subevent - 1) + 5)
            warm_para(subevent, 5) = values(8 * (subevent - 1) + 6)
            warm_para(subevent, 6) = values(8 * (subevent - 1) + 7)
            warm_para(subevent, 7) = values(8 * (subevent - 1) + 2)
        end do
        warm_subevent_count = n
    end subroutine

    ! the warm subevents plus one subevent drawn from the prior, relabelled in time order; E1 must
    ! end up at the hypocentre for two or more subevents (for a 1-subevent warm model, whose centroid
    ! is usually off the epicentre, the new subevent is drawn AT the hypocentre and has to come first).
    !
    ! A crowded warm model gets a second kind of draw. The draw from the prior needs a free slot of
    ! 2 t_gap between two warm subevents, and a model whose pulses follow each other closely has none:
    ! the new subevent can then only be appended after the last one. Chile's 4-subevent model of
    ! 2026-11-17 (centroids at 3.6, 6.1, 9.8 and 14.2 s) was such a model. Its 5-subevent search
    ! started at 0.8716, above the 0.8599 of four subevents, and ended at 0.8714, while one island
    ! of chains of another seed found five subevents at 0.7819, with the new subevent between the
    ! second and the third: the first run selected 3 subevents, the second 5. For a warm model
    ! without a free slot, half of the draws therefore insert the new subevent between two
    ! consecutive warm subevents (interpolated in time, position and depth, plus a jitter) and move
    ! the centroid times of the others apart until every gap is t_gap again. From that 4-subevent
    ! model these draws started at 0.850-0.854 and the 5-subevent searches ended at 0.7795, 0.7801
    ! and 0.7827 (seeds 0, 1000 and 1500, 32 ranks), every island of chains below 0.786 before the
    ! last resampling round. In the full reruns the two searches agree within 0.6% at every N and
    ! both select 4 (0.8124 and 0.8102). A warm model with a free slot is drawn as before, with the
    ! same random numbers: the 3-subevent California search wrote byte-identical chains
    ! (PITFALLS D7).
    subroutine draw_warm_start_model(new_para, ok)
        real, intent(out) :: new_para(max_subevents, 7)
        logical, intent(out) :: ok
        integer :: attempt, n, m, slot, position_draw, left
        real :: swap(7), x_draw, y_draw, inserted(7), fraction, jitter
        logical :: has_epicentral_subevent, warm_model_is_crowded
        has_epicentral_subevent = any(abs(warm_para(1:warm_subevent_count, 1)) < 1e-3 .and. &
                                      abs(warm_para(1:warm_subevent_count, 2)) < 1e-3)
        warm_model_is_crowded = has_epicentral_subevent .and. warm_subevent_count >= 2
        if (warm_model_is_crowded) warm_model_is_crowded = &
            all(warm_para(2:warm_subevent_count, 7) - warm_para(1:warm_subevent_count - 1, 7) < 2 * t_gap)
        ok = .false.
        do attempt = 1, 20000
            new_para = 0
            if (warm_model_is_crowded) then
                if (uniform_deviate(random_state) < 0.5) then
                    left = min(1 + int(uniform_deviate(random_state) * (warm_subevent_count - 1)), warm_subevent_count - 1)
                    fraction = uniform_deviate(random_state)
                    ! jitter: a quarter of the distance a rupture front covers between the two centroid times
                    jitter = 0.25 * vmax_cau * (warm_para(left + 1, 7) - warm_para(left, 7))
                    inserted = warm_para(left, :) + fraction * (warm_para(left + 1, :) - warm_para(left, :))
                    inserted(1) = inserted(1) + (2 * uniform_deviate(random_state) - 1) * jitter
                    inserted(2) = inserted(2) + (2 * uniform_deviate(random_state) - 1) * jitter
                    inserted(3) = max(z_min, min(z_max, inserted(3) + (2 * uniform_deviate(random_state) - 1) * jitter))
                    inserted(4) = dura_min + uniform_deviate(random_state) * (min(inserted(7) * 2, dura_max) - dura_min)
                    inserted(5) = 0.01
                    inserted(6) = 0
                    if (inserted(1) < x_min .or. inserted(1) > x_max .or. inserted(2) < y_min .or. inserted(2) > y_max) cycle
                    if (uniform_deviate(random_state) >= seismicity_prior(inserted(1), inserted(2), grid_spacing, &
                                                                          grid_cell_count, grid_x, grid_y, grid_density)) cycle
                    new_para(1:left, :) = warm_para(1:left, :)
                    new_para(left + 1, :) = inserted
                    new_para(left + 2:neq, :) = warm_para(left + 1:warm_subevent_count, :)
                    do m = left, 1, -1                             ! the earlier subevents move earlier
                        new_para(m, 7) = min(new_para(m, 7), new_para(m + 1, 7) - t_gap)
                    end do
                    do m = left + 2, neq                           ! the later subevents move later
                        new_para(m, 7) = max(new_para(m, 7), new_para(m - 1, 7) + t_gap)
                    end do
                    do m = 1, neq                                  ! no moment before the origin time
                        new_para(m, 4) = min(new_para(m, 4), 2 * new_para(m, 7))
                    end do
                    call place_first_subevent(new_para)
                    if (satisfies_hard_rules(new_para)) then
                        ok = .true.
                        return
                    end if
                    cycle
                end if
            end if
            new_para(1:warm_subevent_count, :) = warm_para(1:warm_subevent_count, :)
            n = neq
            if (has_epicentral_subevent) then
                do position_draw = 1, 10000
                    x_draw = x_min + uniform_deviate(random_state) * (x_max - x_min)
                    y_draw = y_min + uniform_deviate(random_state) * (y_max - y_min)
                    if (uniform_deviate(random_state) < seismicity_prior(x_draw, y_draw, grid_spacing, grid_cell_count, &
                                                                       grid_x, grid_y, grid_density)) exit
                end do
            else
                x_draw = 0
                y_draw = 0
            end if
            new_para(n, 1) = x_draw
            new_para(n, 2) = y_draw
            new_para(n, 3) = z_min + uniform_deviate(random_state) * (z_max - z_min)
            new_para(n, 7) = cen_min + uniform_deviate(random_state) * (cen_max - cen_min)
            new_para(n, 4) = dura_min + uniform_deviate(random_state) * (min(new_para(n, 7) * 2, dura_max) - dura_min)
            new_para(n, 5) = 0.01
            new_para(n, 6) = 0
            do m = 1, neq - 1                                     ! insertion into time order
                do slot = neq, m + 1, -1
                    if (new_para(slot, 7) < new_para(slot - 1, 7)) then
                        swap = new_para(slot, :); new_para(slot, :) = new_para(slot - 1, :); new_para(slot - 1, :) = swap
                    end if
                end do
            end do
            if (neq > 1 .and. (abs(new_para(1, 1)) > 1e-3 .or. abs(new_para(1, 2)) > 1e-3)) cycle
            call place_first_subevent(new_para)           ! the depth of the new E1, or of an E1 of an older run
            if (satisfies_hard_rules(new_para)) then
                ok = .true.
                return
            end if
        end do
    end subroutine

    ! E1 of a model with two or more subevents sits at the hypocentre: at the epicentre and, when the
    ! data headers carry it, at the hypocentre depth
    subroutine place_first_subevent(p)
        real, intent(inout) :: p(max_subevents, 7)
        if (neq <= 1) return
        p(1, 1:2) = 0
        if (first_depth_is_pinned) p(1, 3) = hypocentre_depth
    end subroutine

    ! the time at which the rupture front leaves the hypocentre: the start of E1's moment rate, half
    ! a duration before its centroid
    real function rupture_start_time(p)
        real, intent(in) :: p(max_subevents, 7)
        rupture_start_time = p(1, 7) - p(1, 4) / 2
    end function

    ! the hard rules of the prior; the density factor is applied in the acceptance test
    logical function satisfies_hard_rules(p)
        real, intent(in) :: p(max_subevents, 7)
        real :: v1(2), v2(2), v3(2), angle_continuing, angle_bilateral
        integer :: n
        satisfies_hard_rules = .false.
        do n = 1, neq
            if (p(n, 1) < x_min .or. p(n, 1) > x_max .or. p(n, 2) < y_min .or. p(n, 2) > y_max) return
            if (p(n, 3) < z_min .or. p(n, 3) > z_max .or. p(n, 4) < dura_min .or. p(n, 4) > dura_max) return
            if (p(n, 7) < cen_min .or. p(n, 7) > cen_max) return
            if (p(n, 7) < p(n, 4) / 2) return                                      ! no moment before the origin
            if (n < neq) then
                if (p(n + 1, 7) - p(n, 7) < t_gap - 1e-5) return                     ! labels in time order
            end if
            if (n >= 2) then                                                       ! rupture causality
                if (p(n, 7) <= p(1, 7)) return
                ! no centroid ahead of a front that leaves E1 with vmax_cau when E1 starts
                if (sqrt(sum((p(n, 1:3) - p(1, 1:3))**2)) > vmax_cau * (p(n, 7) - rupture_start_time(p))) return
            end if
            if (n >= 3) then
                ! each new subevent continues the rupture direction from E(n-1) or lies on the far
                ! side of E1 (bilateral): both angles above 45 deg is a zig-zag and is forbidden
                v1 = p(n, 1:2) - p(n - 1, 1:2)
                v2 = p(n, 1:2) - p(1, 1:2)
                v3 = p(n - 1, 1:2) - p(1, 1:2)
                angle_continuing = 90
                angle_bilateral = 90
                if (norm2(v1) > 0 .and. norm2(v3) > 0) &
                    angle_continuing = acos(max(-1.0, min(1.0, dot_product(v1, v3) / (norm2(v1) * norm2(v3))))) * degrees_per_radian
                if (norm2(v2) > 0 .and. norm2(v3) > 0) &
                    angle_bilateral = acos(max(-1.0, min(1.0, -dot_product(v3, v2) / (norm2(v2) * norm2(v3))))) * degrees_per_radian
                if (angle_continuing > 45 .and. angle_bilateral > 45) return
            end if
        end do
        satisfies_hard_rules = .true.
    end function

    ! one proposal into para_proposed. move_kind: 1 x, 2 y, 3 depth, 4 duration, 7 centroid time,
    ! 8 adaptive joint move, 9 adaptive one-dimensional move along a learned direction,
    ! 10 relocation of one subevent to a position drawn from the prior.
    ! Why relocations (5% of proposals): local moves of a few km cannot carry a subevent between
    ! modes tens of km apart - in the published California 3-subevent run chains disagreed on E3's x
    ! (split R-hat 11), and in the 2026 Panama clean-environment run 1 of 16 chains found the best
    ! 2-subevent mode, the median chain stayed 17% worse, and the selected N rested on that one chain.
    subroutine propose(move_kind_out, subevent_out)
        integer, intent(out) :: move_kind_out, subevent_out
        integer :: m, n, first_movable
        para_proposed = para
        subevent_out = 0
        u = uniform_deviate(random_state)
        move_kind_out = 0
        if (u < 0.05) then
            n = 1
            if (neq > 1) n = 2 + min(neq - 2, int(uniform_deviate(random_state) * (neq - 1)))
            move_kind_out = 10
            subevent_out = n
            call draw_relocation(n)
            return
        end if
        u = (u - 0.05) / 0.95
        if (adaptive_ready == 1) then
            if (u < 0.05) then
                move_kind_out = 8
            else if (u < 0.25) then
                move_kind_out = 9
            else
                u = (u - 0.25) / 0.75
            end if
        end if
        if (move_kind_out == 8 .or. move_kind_out == 9) then
            do m = 1, adaptive_count
                gaussian_draws(m) = gaussian_deviate(random_state)
            end do
            if (move_kind_out == 9) step_length = gaussian_deviate(random_state) * exp(0.5 * direction_log_scale)
            do m = 1, adaptive_count
                direction(m) = dot_product(adaptive_factor(m, 1:m), gaussian_draws(1:m))
                if (move_kind_out == 8) then
                    para_proposed(adaptive_subevent(m), adaptive_column(m)) = para(adaptive_subevent(m), adaptive_column(m)) &
                                                                              + direction(m) * exp(0.5 * joint_log_scale)
                else
                    para_proposed(adaptive_subevent(m), adaptive_column(m)) = para(adaptive_subevent(m), adaptive_column(m)) &
                                                                              + step_length * direction(m)
                end if
            end do
            return
        end if
        ! single-parameter move of a uniformly chosen subevent; x, y and (when the data carry the
        ! hypocentre depth) the depth of E1 stay at the hypocentre for two or more subevents
        first_movable = 1
        if (neq > 1 .and. (u < 0.4 .or. (u < 0.6 .and. first_depth_is_pinned))) first_movable = 2
        n = first_movable + min(neq - first_movable, int(uniform_deviate(random_state) * (neq - first_movable + 1)))
        subevent_out = n
        if (u < 0.2) then
            move_kind_out = 1
            para_proposed(n, 1) = para(n, 1) + gaussian_deviate(random_state) * width_x
        else if (u < 0.4) then
            move_kind_out = 2
            para_proposed(n, 2) = para(n, 2) + gaussian_deviate(random_state) * width_y
        else if (u < 0.6) then
            move_kind_out = 3
            para_proposed(n, 3) = para(n, 3) + gaussian_deviate(random_state) * width_z
        else if (u < 0.8) then
            move_kind_out = 4
            para_proposed(n, 4) = para(n, 4) + gaussian_deviate(random_state) * width_duration
        else
            move_kind_out = 7
            para_proposed(n, 7) = para(n, 7) + gaussian_deviate(random_state) * width_centroid
        end if
    end subroutine

    ! relocation of subevent n into para_proposed: position uniform in the search box, kept with
    ! probability equal to the density (normalized to a maximum of 1), depth uniform; time and
    ! duration unchanged
    subroutine draw_relocation(n)
        integer, intent(in) :: n
        integer :: draw_attempt
        real :: x_draw, y_draw
        do draw_attempt = 1, 10000
            x_draw = x_min + uniform_deviate(random_state) * (x_max - x_min)
            y_draw = y_min + uniform_deviate(random_state) * (y_max - y_min)
            if (uniform_deviate(random_state) < seismicity_prior(x_draw, y_draw, grid_spacing, grid_cell_count, &
                                                               grid_x, grid_y, grid_density)) exit
        end do
        para_proposed(n, 1) = x_draw
        para_proposed(n, 2) = y_draw
        para_proposed(n, 3) = z_min + uniform_deviate(random_state) * (z_max - z_min)
    end subroutine

    ! Robbins-Monro adaptation of a log proposal scale toward a target acceptance rate
    subroutine adapt_scale(log_scale, proposal_count, acceptance, target)
        real, intent(inout) :: log_scale
        integer, intent(in) :: proposal_count
        real, intent(in) :: acceptance, target
        log_scale = log_scale + (acceptance - target) / (20.0 + proposal_count)**0.66
    end subroutine

    ! running mean and covariance of the chain (Welford), refactored every 200 states
    subroutine update_adaptive_statistics()
        integer :: m, n
        integer :: factor_ok
        do m = 1, adaptive_count
            adaptive_vector(m) = para(adaptive_subevent(m), adaptive_column(m))
        end do
        adaptive_state_count = adaptive_state_count + 1
        adaptive_delta(1:adaptive_count) = adaptive_vector(1:adaptive_count) - adaptive_mean(1:adaptive_count)
        adaptive_mean(1:adaptive_count) = adaptive_mean(1:adaptive_count) + adaptive_delta(1:adaptive_count) / adaptive_state_count
        do m = 1, adaptive_count
            do n = 1, adaptive_count
                adaptive_m2(m, n) = adaptive_m2(m, n) + adaptive_delta(m) * (adaptive_vector(n) - adaptive_mean(n))
            end do
        end do
        if (adaptive_state_count >= adaptive_start .and. mod(int(adaptive_state_count), 200) == 0) then
            ! Haario scale 2.38^2/d, plus a small floor so a collapsed direction stays proposable
            do m = 1, adaptive_count
                do n = 1, adaptive_count
                    adaptive_factor(m, n) = (2.38**2 / adaptive_count) * adaptive_m2(m, n) / max(adaptive_state_count, 1.0)
                end do
                adaptive_factor(m, m) = adaptive_factor(m, m) + (0.1 * adaptive_width(m))**2
            end do
            call cholesky_lower(adaptive_factor, adaptive_count, max_adaptive, factor_ok)
            if (factor_ok == 1) adaptive_ready = 1
        end if
    end subroutine

    ! Resampling of the population before a burn-in stage at temperature t_new: every chain draws a
    ! parent with probability ~ exp(-(phi_j - phi_min) / (2 sigma^2 phi_ref t_new)), the density of
    ! its current state under the stage's target, and copies the parent's state, adaptive statistics
    ! included. At t_new = 1 a chain 2% worse than the best keeps ~e^-4 of the copies of the best one,
    ! as the posterior itself would; the v2 consolidation kept every chain within 10% with equal
    ! weight, which is why its ensembles mixed separate modes (California 3sub E3 x: split R-hat 9.7).
    ! Collective: every rank calls it at the same iteration.
    subroutine resample_population(t_old, t_new, is_global_round)
        real, intent(in) :: t_old, t_new
        logical, intent(in) :: is_global_round
        real :: weights(rank_count), population_min, cumulative, draw
        logical :: is_candidate(rank_count)
        integer :: parent, member, m
        if (t_old < t_new) stop 'finv: resampling expects a non-increasing temperature'
        if (is_global_round) call write_island_bests()
        call pack_state(own_state)
        call MPI_ALLGATHER(own_state, state_size, MPI_REAL, population, state_size, MPI_REAL, MPI_COMM_WORLD, mpi_status_code)
        ! Islands keep several modes alive through the burn-in: a single population collapsed onto
        ! whichever mode led early (California 3sub, seeds 0/500: best 1.871 vs 1.933, every chain
        ! within 1-2% of its own population's best). The last round is global. Weights are measured
        ! from the best of the candidates, so they cannot all underflow to zero.
        do member = 1, rank_count
            is_candidate(member) = is_global_round .or. (member - 1) / island_size == rank / island_size
        end do
        population_min = minval(population(1, :), mask=is_candidate)
        do member = 1, rank_count
            weights(member) = 0
            if (is_candidate(member)) weights(member) = exp(-(population(1, member) - population_min) &
                                                            / (2 * sigma_relative**2 * like_reference * t_new))
        end do
        draw = uniform_deviate(random_state) * sum(weights)
        cumulative = 0
        parent = rank + 1
        do member = 1, rank_count
            cumulative = cumulative + weights(member)
            if (weights(member) > 0 .and. draw <= cumulative) then
                parent = member
                exit
            end if
        end do
        call unpack_state(population(:, parent))
        if (like < like_best) then
            like_best = like
            para_best = para
        end if
        like_reference = min(like_reference, minval(population(1, :)))
        if (rank == 0) then
            m = count(weights > 0.01 * maxval(weights))
            write(*, *) 'population resampling: best', population_min, ' chains with weight > 1% of max:', m, 'of', &
                        rank_count, merge(' (global)', ' (island)', is_global_round)
        end if
    end subroutine

    ! island_best_misfits.txt (one line per island: island, chains, best misfit so far), written just before
    ! the global round. Until that round the islands are independent populations, so their bests are
    ! independent searches with island_size chains each; after it every chain descends from the best
    ! island and the chains no longer measure how reproducible the search is (their bests then agree
    ! within 1-2% even when two seeds of the same run differ by 5%: California 3sub, seeds 0 and 500,
    ! best 1.761 and 1.847). plot_lcurve.py lists these bests next to the L-curve.
    ! Collective: called by every rank inside resample_population.
    subroutine write_island_bests()
        integer, parameter :: island_unit = 23
        real :: chain_bests(rank_count), own_best(1)
        integer :: island, first_member, last_member
        own_best(1) = like_best          ! an array: gfortran rejects a scalar where other calls pass arrays
        call MPI_ALLGATHER(own_best, 1, MPI_REAL, chain_bests, 1, MPI_REAL, MPI_COMM_WORLD, mpi_status_code)
        if (rank /= 0) return
        open(island_unit, file='island_best_misfits.txt', status='replace')
        do island = 1, (rank_count + island_size - 1) / island_size
            first_member = (island - 1) * island_size + 1
            last_member = min(island * island_size, rank_count)
            write(island_unit, '(I4, I4, F12.5)') island, last_member - first_member + 1, &
                                                  minval(chain_bests(first_member:last_member))
        end do
        close(island_unit)
    end subroutine

    subroutine pack_state(state)
        real, intent(out) :: state(state_size)
        integer :: position
        state = 0
        state(1) = like
        state(2:1 + 7 * max_subevents) = reshape(para, (/ 7 * max_subevents /))
        position = 2 + 7 * max_subevents
        state(position) = adaptive_state_count
        state(position + 1) = real(adaptive_ready)
        state(position + 2) = joint_log_scale
        state(position + 3) = direction_log_scale
        position = position + 4
        state(position:position + max_adaptive - 1) = adaptive_mean
        position = position + max_adaptive
        state(position:position + max_adaptive - 1) = adaptive_width
        position = position + max_adaptive
        state(position:position + max_adaptive**2 - 1) = reshape(adaptive_m2, (/ max_adaptive**2 /))
        position = position + max_adaptive**2
        state(position:position + max_adaptive**2 - 1) = reshape(adaptive_factor, (/ max_adaptive**2 /))
    end subroutine

    subroutine unpack_state(state)
        real, intent(in) :: state(state_size)
        integer :: position
        like = state(1)
        para = reshape(state(2:1 + 7 * max_subevents), (/ max_subevents, 7 /))
        position = 2 + 7 * max_subevents
        adaptive_state_count = state(position)
        adaptive_ready = nint(state(position + 1))
        joint_log_scale = state(position + 2)
        direction_log_scale = state(position + 3)
        position = position + 4
        adaptive_mean = state(position:position + max_adaptive - 1)
        position = position + 2 * max_adaptive               ! widths are identical on every rank
        adaptive_m2 = reshape(state(position:position + max_adaptive**2 - 1), (/ max_adaptive, max_adaptive /))
        position = position + max_adaptive**2
        adaptive_factor = reshape(state(position:position + max_adaptive**2 - 1), (/ max_adaptive, max_adaptive /))
    end subroutine

    subroutine print_progress()
        integer :: n, m_index
        write(*, *) 'sample', iteration, '/', burn_in + nsample, rank
        write(*, *) 'number of sub-events:', neq
        write(*, *) 'like', like, 'line_min', like_best, 'temperature', temperature
        do n = 1, neq
            write(*, '(7F9.3)') para_best(n, 7), para_best(n, 1), para_best(n, 2), para_best(n, 4), para_best(n, 5), &
                                para_best(n, 6), para_best(n, 3)
        end do
        if (iteration > burn_in) then
            write(*, *) 'Accept rate', (100.0 * single_accepted(m_index) / max(1, single_proposed(m_index)), m_index = 1, 7)
            write(*, *) 'AM joint-move accept rate', 100.0 * joint_accepted / max(joint_proposed, 1), ' of ', joint_proposed, &
                        ' scale ', exp(0.5 * joint_log_scale)
            write(*, *) 'AM 1-D dir accept rate', 100.0 * direction_accepted / max(direction_proposed, 1), ' of ', &
                        direction_proposed
            write(*, *) 'relocation accept rate', 100.0 * relocation_accepted / max(relocation_proposed, 1), ' of ', &
                        relocation_proposed
            write(*, *) 'proposals rejected by a hard rule (%)', 100.0 * infeasible_count / max(iteration, 1)
        end if
    end subroutine

    ! one line per sample: neq, ' <neq ', then misfit cen x y dur vr theta depth for every subevent
    ! (format read by plot_lcurve.py, posterior_summary.py, step3/step4 and calibrate_from_1sub.sh)
    subroutine write_sample_line(unit_number, sample_index)
        integer, intent(in) :: unit_number, sample_index
        integer :: n
        write(unit_number, '(I2,A8)', advance='no') neq, ' <neq '
        do n = 1, neq
            write(unit_number, '(8F10.4)', advance='no') samples(sample_index, n, :)
        end do
        write(unit_number, *)
    end subroutine

    ! IEEE-754 test that -Ofast cannot fold away: exponent bits all set = NaN or Inf
    logical function is_finite_real(value)
        real, intent(in) :: value
        is_finite_real = iand(ishft(transfer(value, 0), -23), 255) /= 255
    end function

end program finv

! model vector for forward_model.c: per subevent centroid time, x, y, duration, vr, theta, depth
subroutine model_to_vector(neq, para, model_vector)
    implicit none
    integer, intent(in) :: neq
    real, intent(in) :: para(9, 7)
    real, intent(out) :: model_vector(63)
    integer :: i
    model_vector = 0
    do i = 1, neq
        model_vector((i - 1) * 7 + 1:(i - 1) * 7 + 7) = (/ para(i, 7), para(i, 1), para(i, 2), para(i, 4), para(i, 5), &
                                                           para(i, 6), para(i, 3) /)
    end do
end subroutine

! lower Cholesky factor in place (upper triangle zeroed); ok = 0 when not positive definite
subroutine cholesky_lower(a, n, leading_dimension, ok)
    implicit none
    integer, intent(in) :: n, leading_dimension
    real, intent(inout) :: a(leading_dimension, leading_dimension)
    integer, intent(out) :: ok
    integer :: i, j, k
    real :: s
    ok = 1
    do j = 1, n
        s = a(j, j) - sum(a(j, 1:j - 1)**2)
        if (s <= 0) then
            ok = 0
            return
        end if
        a(j, j) = sqrt(s)
        do i = j + 1, n
            s = a(i, j)
            do k = 1, j - 1
                s = s - a(i, k) * a(j, k)
            end do
            a(i, j) = s / a(j, j)
        end do
    end do
    do j = 1, n
        a(1:j - 1, j) = 0
    end do
end subroutine

! aftershock density at (x, y): bilinear interpolation on the regular grid of seisdens.dat (x
! varying fastest, spacing dx), zero outside the grid
real function seismicity_prior(x, y, dx, cell_count, grid_x, grid_y, grid_density)
    implicit none
    integer, intent(in) :: cell_count
    real, intent(in) :: x, y, dx, grid_x(cell_count), grid_y(cell_count), grid_density(cell_count)
    integer :: columns, rows, ix, iy
    real :: fx, fy, corner(4)
    columns = nint((grid_x(cell_count) - grid_x(1)) / dx) + 1
    rows = cell_count / columns
    ix = floor((x - grid_x(1)) / dx)
    iy = floor((y - grid_y(1)) / dx)
    fx = (x - (grid_x(1) + ix * dx)) / dx
    fy = (y - (grid_y(1) + iy * dx)) / dx
    corner = 0
    if (inside(ix, iy)) corner(1) = grid_density(iy * columns + ix + 1)
    if (inside(ix + 1, iy)) corner(2) = grid_density(iy * columns + ix + 2)
    if (inside(ix, iy + 1)) corner(3) = grid_density((iy + 1) * columns + ix + 1)
    if (inside(ix + 1, iy + 1)) corner(4) = grid_density((iy + 1) * columns + ix + 2)
    seismicity_prior = corner(1) * (1 - fx) * (1 - fy) + corner(2) * fx * (1 - fy) + corner(3) * (1 - fx) * fy + corner(4) * fx * fy
contains
    logical function inside(column, row)
        integer, intent(in) :: column, row
        inside = column >= 0 .and. column <= columns - 1 .and. row >= 0 .and. row <= rows - 1
    end function
end function

! Gaussian deviate (Numerical Recipes polar Box-Muller; the second deviate is discarded as before)
real function gaussian_deviate(state)
    implicit none
    integer, intent(inout) :: state
    real :: v1, v2, r
    real, external :: uniform_deviate
    do
        v1 = 2 * uniform_deviate(state) - 1
        v2 = 2 * uniform_deviate(state) - 1
        r = v1**2 + v2**2
        if (r < 1 .and. r > 0) exit
    end do
    gaussian_deviate = v2 * sqrt(-2 * log(r) / r)
end function

! Numerical Recipes ran3 (Knuth subtractive generator); seeded on the first call from state
real function uniform_deviate(state)
    implicit none
    integer, intent(inout) :: state
    integer, parameter :: mbig = 1000000000, mseed = 161803478, mz = 0
    real, parameter :: fac = 1. / mbig
    integer, save :: iff = 0, inext, inextp, ma(55)
    integer :: i, ii, k, mj, mk
    if (state < 0 .or. iff == 0) then
        iff = 1
        mj = mod(mseed - iabs(state), mbig)
        ma(55) = mj
        mk = 1
        do i = 1, 54
            ii = mod(21 * i, 55)
            ma(ii) = mk
            mk = mj - mk
            if (mk < mz) mk = mk + mbig
            mj = ma(ii)
        end do
        do k = 1, 4
            do i = 1, 55
                ma(i) = ma(i) - ma(1 + mod(i + 30, 55))
                if (ma(i) < mz) ma(i) = ma(i) + mbig
            end do
        end do
        inext = 0
        inextp = 31
        state = 1
    end if
    inext = inext + 1
    if (inext == 56) inext = 1
    inextp = inextp + 1
    if (inextp == 56) inextp = 1
    mj = ma(inext) - ma(inextp)
    if (mj < mz) mj = mj + mbig
    ma(inext) = mj
    uniform_deviate = mj * fac
end function
