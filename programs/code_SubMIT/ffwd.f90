! ffwd.f90 - forward model of one subevent model: waveforms, moment tensors and per-station residuals.
!
! Reads Input.model (one line per subevent: centroid time, x east, y north, duration, rupture
! velocity, rupture direction, depth), evaluates it with the same forward code that finv samples
! (forward_model.c) and lets that code write fm.dat, sta_residual.dat, stf/*.sac and waveforms/*.sac
! for the figures (step 4), the model-based station screen and the weight balancing (step 3 hooks).
!
! Until v2 this program kept its own copy of the forward code and of every data array; the copies
! drifted apart (array sizes that differed between them once overwrote the SH filter
! coefficients). It now holds nothing but the model.
!
! Model lists (2026-11-16). Given a file name, ffwd evaluates every model in that file (models
! separated by blank lines, each in the Input.model format) after ONE read of the data and the
! Green's functions. The read takes about 20 s on the 2026 test machine and one model well under a
! second, so evaluate_misfit_components.sh scores the best sample of all 32 chains in one run.
!
! Run in an inversion or forward directory after "make ffwd":   ./ffwd
! FFWD_COMPACT=1 ./ffwd prints the misfit lines only (used by the weight balancing).
! FFWD_COMPACT=1 ./ffwd MODELS_FILE prints them for every model of MODELS_FILE.
program ffwd
    implicit none
    integer, parameter :: max_subevents = 9
    real :: model(7, max_subevents), penalized_misfit, source_s_velocity
    integer :: subevent_count, read_status, evaluated_model_count
    character(len=512) :: line, model_file_name
    logical :: is_model_list

    ! Input.model holds one model, which ends at its first blank line; a model list holds several
    is_model_list = command_argument_count() >= 1
    model_file_name = 'Input.model'
    if (is_model_list) call get_command_argument(1, model_file_name)
    open(10, file=trim(model_file_name), status='old')

    call sub_init(source_s_velocity)
    call enable_forward_outputs()
    evaluated_model_count = 0
    subevent_count = 0
    do
        read(10, '(A)', iostat=read_status) line
        if (read_status == 0 .and. len_trim(line) > 0) then
            if (subevent_count == max_subevents) error stop 'ffwd: a model has more than 9 subevents'
            subevent_count = subevent_count + 1
            read(line, *) model(:, subevent_count)
            cycle
        end if
        ! a blank line or the end of the file closes the model read so far
        if (subevent_count > 0) then
            call sub_forward(subevent_count, model, penalized_misfit)
            write(*, '(A, F12.6)') 'penalized misfit', penalized_misfit
            evaluated_model_count = evaluated_model_count + 1
            subevent_count = 0
        end if
        if (read_status /= 0 .or. .not. is_model_list) exit
    end do
    close(10)
    if (evaluated_model_count == 0) error stop 'ffwd: the model file holds no model'
end program ffwd
