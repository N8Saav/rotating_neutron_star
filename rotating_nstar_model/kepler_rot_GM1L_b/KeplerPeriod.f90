program FreqToPeriod
    implicit none
    real :: a, b, c, freq, period, pi = ACOS(-1.0), freqhz
    integer :: n

    open(100, file="kepler_frequency_GM1L_b.dat", status='old')
    open(200, file="kepler_period_GM1L_b.dat", status='unknown')

    read(100,2)
    2 format (/) 

    do n =1, 40
        read(100,*) a, b, c, freq, a, a, a, a, a
        freqhz = freq/(2*pi)
        period = (1. / freqhz) * (1.0e03)
        write(200,*) n, period
        write(*,*) period
    end do

    close(100)
    close(200)
    
end program FreqToPeriod