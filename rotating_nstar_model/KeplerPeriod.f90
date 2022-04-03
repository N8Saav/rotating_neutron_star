program FreqToPeriod
    implicit none
    real :: a, b, c, freq, period
    integer :: n

    open(100, file="GM1L_h_KeplerFreq.dat", status='old')
    open(200, file="GM1L_h_KeplerPeriod.dat", status='unknown')

    read(100,2)
    2 format (/) 

    do n =1, 41
        read(100,*) a, b, c, freq, a, a, a, a, a
        period = (1. / freq) * (10**(3))
        write(200,*) n, period
        write(*,*) period
    end do

    close(100)
    close(200)
    
end program FreqToPeriod