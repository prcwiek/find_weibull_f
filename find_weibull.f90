! Simple program for finding the Weibull distribution parameters
! k shape factor and c scale factor.
! The file with the example wind measurement data set comes
! from the measurement mast US Virgin Islands St. Thomas Bovoni and
! was downloaded from the site 
! https://midcdmz.nrel.gov/apps/sitehome.pl?site=USVILONA.
! 
! Roberts, O.; Andreas, A.; (1997). United States Virgin Islands:
! St. Thomas & St. Croix (Data); NREL Report No. DA-5500-64451.
! http://dx.doi.org/10.7799/1183464 
! https://midcdmz.nrel.gov/
!
! https://midcdmz.nrel.gov/
!
! Sorting function from 
! https://fortran-lang.discourse.group/t/modern-fortran-sample-code/2019/4
!

program find_weibull
    
    implicit none
        
    integer, parameter :: maxread = 10**6

    logical :: file_exists
    
    character(len=9999) :: fname

    integer :: i, nread
    integer :: funit, ierr
    integer :: niter
    
    real :: k ! shape factor
    real :: c ! scale factor
    real :: ws_mean, ws_median
    real :: kmax, kmin
    real :: eps
    real, allocatable :: rtemp(:), ws(:)

    ! check if an argument is present
    if(command_argument_count() < 1) stop 'Usage: find_weibull.out <filename>'
    
    ! get a filename from an argument
    call get_command_argument(1, fname)

    ! check if file exists
    inquire(file = fname, exist = file_exists)
    if(.not. file_exists) stop 'File does not exist!'
 
    allocate(rtemp(maxread))
    nread = maxread
    open(newunit = funit, file = fname, action = "read", status = "old")
    read(funit, *)      ! read the header
    do i = 1, maxread
            read(funit, *, iostat = ierr) rtemp(i)
            if (ierr /= 0) then ! reached end of file
                    nread = i - 1
                    exit
            end if
    end do
    close(funit)

    allocate(ws(nread), source=rtemp(:nread))
    deallocate(rtemp)
    
    ! range for searching k and c
    kmin = 1.0
    kmax = 8.0
    ! accuracy
    eps = 0.000001
    ! number of iterations
    niter = 50
    ! shape factor
    k = bisection(ws, kmin, kmax, eps, niter)
    
    ! mean wind speed 
    ws_mean = sum(ws) / nread
    
    ! scale factor
    c = ws_mean * (0.586 + 0.433/k)**(-1/k)
    
    ! median wind speed
    ws_median = median(ws)

    print*, "Found Weibull distribution parameters:"
    print*
    write(*, "(A,F4.2,A)") " shape factor k: ", k
    write(*, "(A,F4.2,A)") " scale factor c: ", c
    print*
    write(*, "(A,F4.2,A)") " Mean wind speed: ", ws_mean, " m/s"
    write(*, "(A,F4.2,A)") " Median wind speed: ", ws_median, " m/s"

    deallocate(ws)

contains

    ! k estimator
    pure real function k_estimator(x, kin) result(res)
        
        real, intent(in) :: x(:)
        real, intent(in) :: kin
        
        integer :: n
        
        n = size(x)

        res = (sum(x**3) / n) / ((sum(x) / n)**3)
        res = res * (gamma(1.0+1.0/kin)**3) - gamma(1.0 + 3.0/kin)

    end function k_estimator


    real function bisection(x, ikmin, ikmax, eps, iter) result (k)
        
        real, intent(in) :: x(:)
        real, intent(in) :: ikmin, ikmax, eps
        integer, intent(in) :: iter
        !real :: k
        
        integer :: j
        
        real :: fk, fkk
        real :: fkmin, fkmax
        real :: kmin, kmax
        
        ! initial values
        kmin = ikmin
        kmax = ikmax
        fkmin = k_estimator(x, kmin)
        fkmax = k_estimator(x, kmax)
        
        if (fkmin * fkmax > 0) then
                print *, "Error: Both estimated k values are greater than zero!"
                k = 0
                stop
        end if
        
        do j = 1, iter
                k = (kmin + kmax) / 2
                fk = k_estimator(x, k)
                fkk = (fkmax - fkmin) / (kmax - kmin)
                
                
                if (abs(fk/fkk) - eps > 0) then
                        if (fk*fkmin < 0) then
                                kmax = k
                                fkmax = fk
                        else
                                if (fk * fkmin == 0)  then
                                        return
                                end if
                                kmin = k
                                fkmin = fk
                        end if
                else
                        return
                end if
                        
        end do
        k = 0
    end function bisection

    pure real function median(x) result(res)
        
        real, intent(in) :: x(:)
        real :: xs(size(x))
        
        integer :: n
               
        xs = qsort(x)
        
        n = size(x)

        if (mod(n, 2) == 0) then
                res = (xs(n/2) + xs(n/2+1)) / 2
        else
                res = xs(n+1)/2
        end if
    end function median

    ! sorting function from 
    ! https://fortran-lang.discourse.group/t/modern-fortran-sample-code/2019/4
    pure recursive function qsort(data) result(sorted)
        real, intent(in) :: data(:)
        real             :: sorted(size(data))
        
        if (size(data) > 1) then
           sorted = [qsort(pack(data(2:),data(2:)<data(1))), data(1), &
                     qsort(pack(data(2:),data(2:)>=data(1)))]
        else
           sorted = data
        end if
    end function qsort

end program find_weibull
