!simulation of slab with kx=1.5, evanescent waves, truncated at z=0
!first interface at z=1, second at z=3 
!(i.e. distance parametised by d_s, distance from source to interface)


implicit none

double precision :: kxtilde, ztildestepfrac, xtildestepfrac, dsourcetilde, smallval, eta
double precision :: secondinterface, xtildei, ztildei, c=3D8, xtildef, ztildef, PI, eps1, mu1, mu2
complex*16 :: A, C2, D, T, n1, n2, eps2, kz1tilde, kz2tilde, chi, i, denominator
integer*4 :: m, n, p, xtildesize, ztildesize
character :: filename*150, ti*10

double precision, dimension(:), allocatable :: xtildearray
double precision, dimension(:), allocatable :: ztildearray

complex*16, dimension(:,:), allocatable :: Eyrealarray
complex*16, dimension(:,:), allocatable :: Eyarray


!The xtilde values etc. are the real values of x, turned into dimensionless parameters via d_s
! i.e. xtilde = x/d_s etc.
i = (0.0,1.0)
ztildei=0
xtildei=-10
ztildef=20
xtildef=10
ztildestepfrac=0.1
xtildestepfrac=0.01
xtildesize = anint(((xtildef-xtildei)/xtildestepfrac))
ztildesize = anint(((ztildef-ztildei)/ztildestepfrac))

allocate(Eyarray(0:xtildesize,0:ztildesize))
allocate(xtildearray(0:xtildesize))
allocate(ztildearray(0:ztildesize))


!The distance from the source to the interface parametised by d (in zspace)
dsourcetilde = 1

!position of second interface
secondinterface = 3*dsourcetilde
eps1=1.0
mu1=1.0
eps2=-1.0
mu2=-1.0
!ti = '0' !also change thetai
PI=4.D0*DATAN(1.D0) 
!ensures maximum precision on any architechture apparently
!smallval = 0.1
!thetai= (smallval/180.0)*PI !also change ti

n1=SQRT(eps1*mu1)
n2=SQRT(eps2*mu2)


if ((RealPart(n2) < 0 .and. RealPart(eps2) > 0) .or. &
 (RealPart(eps2) < 0 .and. mu2<0 .and. RealPart(n2) > 0 )) then
	n2 = -1 * n2
end if

do p=0, xtildesize
	xtildearray(p)= xtildei + p*xtildestepfrac
end do

do p=0, ztildesize
	ztildearray(p)= ztildei + p*ztildestepfrac
end do
eta=1.0
kxtilde=1.5
kz1tilde=sqrt(abs((n1*eta)**2 - kxtilde**2))
print*, " n2= ", n2 
kz2tilde=sqrt(abs((n2*eta)**2 - kxtilde**2))
chi = (kz2tilde*mu1)/(kz1tilde*mu2)
!Define coefficients

denominator=(((chi+1)**2)-((chi-1)**2)*exp(-4*kz2tilde*dsourcetilde)) 
D=((2*(chi-1)*exp(-kz1tilde*dsourcetilde)*exp(-2*kz2tilde)*dsourcetilde))/denominator
T=(4*chi*exp(-kz1tilde*dsourcetilde)*exp(-2*kz2tilde*dsourcetilde))/denominator
C2=(2*exp(-kz1tilde*dsourcetilde)*(chi+1))/denominator
!A=C2 + D*exp(-2*kz2tilde*dsourcetilde) - exp(-kz1tilde*dsourcetilde)
A=0 !should be zero but floating point is causing problems

do m=0, ztildesize
	do n=0, xtildesize
		Eyarray(n,m)=0
		if(ztildearray(m) <= dsourcetilde) then
			!REGION 1
			Eyarray(n,m) = (exp(-kz1tilde*ztildearray(m))+ &
			 A*exp(kz1tilde*(ztildearray(m)-dsourcetilde)) )*exp(i*kxtilde*xtildearray(n))
			!print*, " eyarray= ", abs(Eyarray(n,m))


		elseif ( ztildearray(m) <= (secondinterface) ) then
			!REGION 2
			Eyarray(n,m) = ( C2*exp(-kz2tilde*(ztildearray(m)-dsourcetilde)) +&
			 D*exp(-kz2tilde* (3*dsourcetilde - ztildearray(m)) ) )*exp(i*kxtilde*xtildearray(n))
			!print*, " eyarray= ", Eyarray(n,m)
		else
			!REGION 3
			Eyarray(n,m) = (T*exp(-kz1tilde*(ztildearray(m)- 3*dsourcetilde)))*exp(i*kxtilde*xtildearray(n))
			!print*, " eyarray= ", Eyarray(n,m)
		end if
	end do	
end do


write(filename,20) 'data/evandoubint',kxtilde,'kx',eta,'eta',secondinterface,'secint.dat'
open(unit=2,file= filename)
	


do m=0, ztildesize
	do n=0, xtildesize
		write(2,10) xtildearray(n), ztildearray(m), abs(Eyarray(n,m)), realpart(Eyarray(n,m)) !intensity so square
	end do
end do
10	format(4e15.5,4e15.5,4e15.5,4e15.5)
20 	format(A,f3.1,A,f3.1,A,f3.1,A)


end


