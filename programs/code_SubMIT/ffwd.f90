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
! Run in an inversion or forward directory after "make ffwd":   ./ffwd
! FFWD_COMPACT=1 ./ffwd prints the misfit lines only (used by the weight balancing).
program ffwd
    implicit none
    integer, parameter :: max_subevents = 9
    real :: model(7, max_subevents), penalized_misfit, source_s_velocity
    integer :: subevent_count, read_status
    character(len=512) :: line

    ! one line per subevent; blank lines end the list
    subevent_count = 0
    open(10, file='Input.model', status='old')
    do
        read(10, '(A)', iostat=read_status) line
        if (read_status /= 0 .or. len_trim(line) == 0) exit
        if (subevent_count == max_subevents) stop 'ffwd: Input.model has more than 9 subevents'
        subevent_count = subevent_count + 1
        read(line, *) model(:, subevent_count)
    end do
    close(10)
    if (subevent_count == 0) stop 'ffwd: Input.model is empty'

    call sub_init(source_s_velocity)
    call enable_forward_outputs()
    call sub_forward(subevent_count, model, penalized_misfit)
    write(*, '(A, F12.6)') 'penalized misfit', penalized_misfit
end program ffwd
