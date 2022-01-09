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

program find_weibull
	implicit none
	
	!character(len=10) :: header_name
	character(len=*), parameter :: fname = "WS125.txt"
	integer, parameter :: maxread = 10**6
	
	integer :: i, nread
	integer :: funit, ierr
	integer :: niter
	
	real :: k ! shape factor
	real :: c ! scale factor
	real :: ws_mean, ws_median
	real :: kmax, kmin
	real :: eps
	real, allocatable :: rtemp(:), ws(:), wss(:)

	allocate(rtemp(maxread))
	nread = maxread
	open(newunit = funit, file = fname, action = "read", status = "old")
	read(funit, *)
	do i = 1, maxread
		read(funit, *, iostat = ierr) rtemp(i)
		if (ierr /= 0) then ! reached end of file
			nread = i - 1
			exit
		end if
	end do
	close(funit)

	allocate(ws(nread))
	allocate(wss(nread))
	ws = rtemp(:nread)
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
	! commented as the sorting function slows down the program
	!ws_median = median(ws)

	print*, "Found Weibull distribution parameters:"
	print*
	write(*, "(A,F4.2,A)") " shape factor k: ", k
	write(*, "(A,F4.2,A)") " scale factor c: ", c
	print*
	write(*, "(A,F4.2,A)") " Mean wind speed: ", ws_mean, " m/s"
	!write(*, "(A,F4.2,A)") " Median wind speed: ", ws_median, " m/s"

	deallocate(ws)

contains

! k estimator
function k_estimator(x, kin) result (kout)
	implicit none
	
	real, dimension(1:), intent(in) :: x
	real, intent(in) :: kin
	
	integer :: i, n
	real :: kout
	real :: xmean, xmean3
	real, allocatable :: x3(:)
	
	n = size(x)
	
	xmean = (sum(x) / n)**3

	allocate(x3(n))

	x3 = x**3
	xmean3 = sum(x3) / n
	kout = xmean3 / xmean
	kout = kout * (gamma(1.0+1.0/kin)**3) - gamma(1.0 + 3.0/kin)

	deallocate(x3)
end function k_estimator


function bisection(x, ikmin, ikmax, eps, iter) result (k)
	implicit none
	
	real, dimension(1:), intent(in) :: x
	real, intent(in) :: ikmin
	real, intent(in) :: ikmax
	real, intent(in) :: eps
	real :: k, fk, fkk
	integer, intent(in) :: iter
	
	
	integer :: j
	
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

function median(x) result(m)
	implicit none
	
	real, dimension(1:), intent(in) :: x
	real, allocatable :: xs(:)
	real :: m
	
	integer :: f, n
	
	n = size(x)
	allocate(xs(n))
	
	xs = sort_si(x)
	
	if (mod(n, 2) == 0) then
		m = (xs(n/2) + xs(n/2+1)) / 2
	else
		m = xs(n+1)/2
	end if
	
	deallocate(xs)
end function median

function sort_si(x) result(xs)
! Simpele soring with straight insertion
	implicit none
	
	real, dimension(1:), intent(in) :: x
	real, allocatable :: xs(:)
	
	integer :: i,j,n
	
	real :: a
	
	n = size(x)
	allocate(xs(n))
	
	xs = x
	
	do j = 2, n
		a = xs(j)
		do i = j-1, 1, -1
			if (xs(i) <= a) then
				goto 99
			end if
			xs(i+1) = xs(i)
		end do
		i = 0
99		xs(i+1) = a
	end do
	
end function sort_si

end program find_weibull
