!singleinterface.f90
! By: Alex McMurray 
!
! Derived from rotated frame, but adjusted back into nonrotated space
! Consider only single interface will later expand to double interface (slab)
! Consider only positive refraction (negative refraction has a backward wave)
! This approach will differ from the congaussian code at high angles of incidence
! These results should be compared once the slab is obtained.

program singleinterface
implicit none !doesn't assume types from names, must be declared explicitly

! Declare stuff here - check these are all necessary
double precision :: x = 0.0, z = 0.0, thetai, eta, xtildestepfrac, ztildestepfrac, kxtildestepfrac, dsourcetilde!, kz2tilde
double precision :: xtildei, ztildei, c=3D8, xtildef, ztildef, Eyout, PI, Eyout2, sigmatilde, temp, eps1, eps2, mu1, mu2, n1, n2
complex*16 :: Ey, Ie, Re, Ce, De, Te, kx, kz1, kz2, A, B, C2, D2, F, G2, H, I2, J, K, i, test, kc, integral
complex*16 :: r, t, kz1tilde, kx1prime, kx2prime, xprime, kz2tilde, zprime, kz1prime, kz2prime
integer*4 :: m, n, p, xtildesize, ztildesize,tilen !, SIZE
character :: filename*80, ti*10
double precision, dimension(:), allocatable :: xtildearray
double precision, dimension(:), allocatable :: ztildearray
complex*16, dimension(0:200) :: kxtildearray
complex*16, dimension(:,:), allocatable :: Eyrealarray
complex*16, dimension(:,:), allocatable :: Eykspacearray

!The xtilde values etc. are the real values of x, turned into dimensionless parameters via the 'thickness' d
! i.e. xtilde = x/d etc.

ztildei=0
xtildei=-10
ztildef=20
xtildef=10
ztildestepfrac=0.1
xtildestepfrac=0.1

dsourcetilde = 10
!The distance from the source to the interface parametised by d (in zspace)


eps1=1.0
mu1=1.0
eps2=5.0
mu2=1.0

n1=SQRT(eps1*mu1)
n2=SQRT(eps2*mu2)

xtildesize = int(((xtildef-xtildei)/xtildestepfrac)) + 1
ztildesize = int(((ztildef-ztildei)/ztildestepfrac)) + 1

PI=4.D0*DATAN(1.D0) 
!ensures maximum precision on any architechture apparently

thetai= PI/4.0 
!angle of incidence (i.e. rotation of frame, but we are staying in non-rotated space)

!allocate arrays with desired shape
allocate(xtildearray(0:xtildesize))
allocate(ztildearray(0:ztildesize))
allocate(Eyrealarray(0:xtildesize,0:ztildesize))
allocate(Eykspacearray(0:xtildesize,0:ztildesize))


i = (0.0,1.0)
!Can just use normal functions as modern fortran can determine the type required, 
!specialist csqrt etc. are obsolete

eta = 1 
!go from 0.1 to 5, dimensionless parameter equal to omega*d/c
! where d is the thickness of the slab and lambda is the free space wavelength of the incident light
! w and d and lambda are all replaced by eta


sigmatilde = 1 
!sigma is form of gaussian beamwidth parameter, here taken to be the squared value, dimensions of area
!sigma tilde is sigma/d^2 so it is dimensionless parameter


ti = '0.16pi'
tilen=LEN(TRIM(ti)) !this is just for filename purposes
kxtildestepfrac = ((eta*COS(thetai))/100.0) !no 2pi term?


do p=0, xtildesize
	xtildearray(p)= xtildei + p*xtildestepfrac
end do

do p=0, ztildesize
	ztildearray(p)= ztildei + p*ztildestepfrac
end do
do p=0, 200
	kxtildearray(p)= -(eta*SIN(PI/2.0)) + p*kxtildestepfrac !used to be 0+ 
end do

do m=0, ztildesize
	do n=0, xtildesize
		Eykspacearray(n,m)=0
		
			do p=1, 199
			!avoid singularities
				kz1tilde=SQRT(n1*eta**2 - kxtildearray(p)**2)
				kz2tilde=SQRT(n2*eta**2 - kxtildearray(p)**2)

				xprime=(-SIN(thetai)*ztildearray(m) + COS(thetai)*xtildearray(n))
				zprime=(COS(thetai)*ztildearray(m) + SIN(thetai)*xtildearray(n))

				kx1prime=(COS(thetai)*kxtildearray(p) - SIN(thetai)*kz1tilde)
				kz1prime=(COS(thetai)*kz1tilde + SIN(thetai)*kxtildearray(p))

				kx2prime=(COS(thetai)*kxtildearray(p) - SIN(thetai)*kz2tilde)
				kz2prime=(COS(thetai)*kz2tilde + SIN(thetai)*kxtildearray(p))


				A=exp(i*kz1tilde*dsourcetilde)
				C=exp(i*kz2tilde*dsourcetilde)
				D2=-(kz1tilde/mu1)*exp(i*kz1tilde*dsourcetilde)
				F=-(kz2tilde/mu2)*exp(i*kz2tilde*dsourcetilde)

				t=(D2**2 + A**2)/(D2*F + A*C2)
				r=A*( (C2*(D2**2) - (A*D2*F))/(D2*F + A*C2) )

				if(ztildearray(m) <= dsourcetilde) then
					integral= ((1.0/sqrt(2*PI))*sqrt(sigmatilde)) &
					*(EXP( (-sigmatilde*((kx1prime)**2)/2.0) & !changed rotational terms to the primes instead of just using the direct tildes
					 + (i* kx1prime * xprime  ) + (i*kz1prime*zprime)) &
			 		 *kxtildestepfrac/COS(thetai))

					Eykspacearray(n,m) = Eykspacearray(n,m) + (integral + r*integral)
				else
					integral= ((1.0/sqrt(2*PI))*sqrt(sigmatilde)) &
					*(EXP( (-sigmatilde*((kx2prime)**2)/2.0) & !as above
					 + (i*kx2prime*xprime) + (i*kz2prime*zprime)) &
			 		 *kxtildestepfrac/COS(thetai))

					Eykspacearray(n,m) = Eykspacearray(n,m) + (t*integral)
				end if


			end do 
	end do 
end do

write(filename,20) 'data/singint',ti(1:tilen),'rads',eta,'eta', sigmatilde,'sigmatilde.dat'
open(unit=2,file= filename)
	


do m=0, ztildesize
	do n=0, xtildesize
		write(2,10) xtildearray(n), ztildearray(m), abs(Eykspacearray(n,m)), realpart(Eykspacearray(n,m))
	end do
end do
10	format(4e15.5,4e15.5,4e15.5,4e15.5)
20 	format(A,A,A,f3.1,A,f3.1,A)


end program singleinterface


