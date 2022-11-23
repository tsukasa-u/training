subroutine count(func, t, N, a, result, runtime)
    
    ! use,intrinsic :: iso_fortran_env
    implicit none

    interface
        function h(t, N, a) result(retval)
            implicit none

            integer(4), intent(in) :: N
            real(4),    intent(in) :: t
            real(4),    intent(in) :: a(0:N)
            real(4) :: retval

        end function h
    end interface

    procedure(h) :: func
    integer(4), intent(in) :: N
    real(4),    intent(in) :: t
    real(4),    intent(in) :: a(0:N)
    real(4), intent(out) :: result
    real(4), intent(out) :: runtime

    real(4) :: time_begin_s, time_end_s

    call cpu_time(time_begin_s)
    result = func(t, N, a)
    call cpu_time(time_end_s)

    runtime = time_end_s - time_begin_s

end subroutine count

program main
    implicit none

    integer,parameter :: n = 100000000 + 1
    real(4) :: a(n)
    real(4) :: t = 0.9
    
    real(4) :: results_f(4)
    real(4) :: runtimes_f(4)
    real(4) :: results_g(4)
    real(4) :: runtimes_g(4)

    integer(4) :: times(4) = [100000, 1000000, 10000000, 100000000]

    integer(4) :: m

    call random_number(a)

    do m = 1, 4, 1
        call count(f, t, times(m), a, results_f(m), runtimes_f(m))
        call count(g, t, times(m), a, results_g(m), runtimes_g(m))
    enddo

    
    do m = 1, 4, 1
        print *, times(m), results_f(m), runtimes_f(m)
    enddo

    print *

    do m = 1, 4, 1
        print *, times(m), results_g(m), runtimes_g(m)
    enddo

    contains

        function f(t, N, a) result(retval)
            implicit none
        
            integer(4), intent(in) :: N
            real(4),    intent(in) :: t
            real(4),    intent(in) :: a(0:N)
            real(4) :: retval
        
            integer(4) :: m
        
            retval = 0.0
            do m = 0, N, 1
                retval = retval +  a(m)*t**m
            enddo
        
        end function f
        
        function g(t, N, a) result(retval)
            implicit none
        
            integer(4), intent(in) :: N
            real(4),    intent(in) :: t
            real(4),    intent(in) :: a(0:N)
            real(4) :: retval
        
            integer(4) :: m
        
            retval = 0.0
            do m = N, 0, -1
                retval = retval*t + a(m)
            enddo
        
        end function g

end program main