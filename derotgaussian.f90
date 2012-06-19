!rotgaussian1.f90
! By: Alex mcMurray and Alun Daley
!
! Try to test bug where angle is incorrect via using a rotated frame
!
! My real parts of the fourier transform approach seem to have some small periodic artefacts on them away from the beam, 
!  but they are not present on any of the other plots.

program rotgaussian1
implicit none !doesn't assume types from names, must be declared explicitly

! Declare stuff here
double precision :: x = 0.0, z = 0.0, thetai, eta, xtildestepfrac, ztildestepfrac, kxprimestepfrac
double precision :: xtildei, ztildei, c=3D8, xtildef, ztildef, Eyout, PI, Eyout2, sigmatilde, temp
complex*16 :: Ey, Ie, Re, Ce, De, Te, kx, kz1, kz2, A, B, C2, D2, F, G2, H, I2, J, K, i, test, kc, integral
integer*4 :: m, n, p, xtildesize, ztildesize,tilen !, SIZE
character :: filename*80, ti*10
double precision, dimension(:), allocatable :: xtildearray
double precision, dimension(:), allocatable :: ztildearray
complex*16, dimension(0:200) :: kxprimearray
!double precision, dimension(:), allocatable :: xprimearray
!double precision, dimension(:), allocatable :: zprimearray
complex*16, dimension(:,:), allocatable :: Eyrealarray
complex*16, dimension(:,:), allocatable :: Eykspacearray
double precision, dimension(:,:,:), allocatable :: rotspacearray
!rotspace stores rotated coordinates of each point in real space

!The xtilde values etc. are the real values of x, turned into dimensionless parameters via the 'thickness' d
! i.e. xtilde = x/d etc.

ztildei=-10
xtildei=-10
ztildef=10
xtildef=10
ztildestepfrac=0.1
xtildestepfrac=0.1

xtildesize = int(((xtildef-xtildei)/xtildestepfrac)) + 1
ztildesize = int(((ztildef-ztildei)/ztildestepfrac)) + 1

PI=4.D0*DATAN(1.D0) ! ensures maximum precision on any architechture apparently

thetai= PI/6.0 !used to work out rotated frame

!allocate arrays with desired shape
allocate(xtildearray(0:xtildesize))
allocate(ztildearray(0:ztildesize))
allocate(Eyrealarray(0:xtildesize,0:ztildesize))
allocate(Eykspacearray(0:xtildesize,0:ztildesize))
!allocate(xprimearray(0:xtildesize))
!allocate(zprimearray(0:ztildesize))
allocate(rotspacearray(0:xtildesize,0:ztildesize,1:2))
i = (0.0,1.0)
!Can just use normal functions as modern fortran can determine the type required, 
!specialist csqrt etc. are obsolete

eta = 1 !go from 0.1 to 5, dimensionless parameter equal to omega*d/c
! where d is the thickness of the slab and lambda is the free space wavelength of the incident light
! w and d and lambda are all replaced by eta

!First test only assumes vaccuum therefore no need for epsilon etc.


sigmatilde = 1 !sigma is form of gaussian beamwidth parameter, here taken to be the squared value, dimensions of area
!sigma tilde is sigma/d^2 so it is dimensionless parameter


ti = '0.16pi'
tilen=LEN(TRIM(ti)) !this is just for filename purposes
kxprimestepfrac = ((eta)/100.0)! can be simplified as is only considering normal (in rotated frame)
!and vacuum.
!no 2pi term? - check this


do p=0, xtildesize
	xtildearray(p)= xtildei + p*xtildestepfrac
end do

do p=0, ztildesize
	ztildearray(p)= ztildei + p*ztildestepfrac
end do
do p=0, 200
	kxprimearray(p)= -(eta) + p*kxprimestepfrac !used to be 0+ 
	!no sin due to normal inc. in rot frame? - check this
	!print *, kxprimearray(p)
end do
!calculuate rotated space co-ordinates from parametised coordinates.

do m=0, ztildesize
	do n=0, xtildesize
		!xprimetilde coordinates:
		rotspacearray(n,m,1)= -sin(thetai)*ztildearray(m) + cos(thetai)*xtildearray(n)
		!zprimetilde coordinates:
		rotspacearray(n,m,2)= cos(thetai)*ztildearray(m) + sin(thetai)*xtildearray(n)
	end do
end do		


! do p=0, xtildesize
! 	xprimearray(p)= -sin(thetai)*ztildearray(p) + cos(thetai)*xtildearray(p)
! 	print *, xprimearray(p)
! end do

! do p=0, ztildesize
! 	zprimearray(p)= cos(thetai)*ztildearray(p) + sin(thetai)*xtildearray(p)
! end do

!Obtain E_source directly from real space equation:

do m=0, ztildesize
	do n=0, xtildesize
		Eyrealarray(n,m) = EXP((-rotspacearray(n,m,1)**2)/(2.0*sigmatilde) + i*eta*rotspacearray(n,m,2))
	end do 
end do

! Obtain E_source from fourier transform for comparison
do m=0, ztildesize
	do n=0, xtildesize
		integral=0
		do p=0, 200
			integral = integral + (EXP( (-sigmatilde*(kxprimearray(p)**2)/2.0) &
			 + (i*kxprimearray(p)*rotspacearray(n,m,1)) + (i*sqrt((eta**2) &
			  - (kxprimearray(p)**2))*rotspacearray(n,m,2) ) )*kxprimestepfrac)
		end do 
		Eykspacearray(n,m) = (((1.0/sqrt(2*PI))*sqrt(sigmatilde))*integral)
	end do 
end do

write(filename,20) 'data/rottest',ti(1:tilen),'rads',eta,'eta', sigmatilde,'sigmatilde.dat'
open(unit=2,file= filename)
	


do m=0, ztildesize
	do n=0, xtildesize
		write(2,10) xtildearray(n), ztildearray(m), abs(Eyrealarray(n,m)), abs(Eykspacearray(n,m)),&
		 realpart(Eyrealarray(n,m)), realpart(Eykspacearray(n,m))
	end do
end do
10	format(4e15.5,4e15.5,4e15.5,4e15.5,4e15.5,4e15.5)
20 	format(A,A,A,f3.1,A,f3.1,A)

!real part of transform has artefacts, check this

end program rotgaussian1