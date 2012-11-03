!doubleinterface.f90
!Fixing major bug due to incorrect integral limits"
!therefore change kxtilde for kxtildeprime in integral
!Single interface code using rotated frame approach but translated
!to normal frame.
!Supports command line interface for variables and calls to MATLAB for plotting.
!COPY FOR WORKING ON

!
implicit none !doesn't assume types from names, must be declared explicitly

! Declare stuff here - check these are all necessary
double precision :: thetai, eta, xtildestepfrac, ztildestepfrac, kxtildeprimestepfrac, dsourcetilde
double precision :: xtildei, ztildei, c=3D8, xtildef, ztildef, PI, sigmatilde, eps1, mu1, mu2, etatest
complex*16 :: A, B, C2, D, i, integral, integral2, chi, n1, n2, eps2
complex*16 :: r, t, kztildeprime, kz1tilde,kz2tilde,kxtilde, Ce, De, rnumerator, rdenominator
integer*4 :: m, n, p, xtildesize, ztildesize,tilen
character :: filename*80, ti*10
double precision, dimension(:), allocatable :: xtildearray
double precision, dimension(:), allocatable :: ztildearray
complex*16, dimension(0:200) :: kxtildeprimearray
complex*16, dimension(:,:), allocatable :: Eyrealarray
complex*16, dimension(:,:), allocatable :: Eykspacearray





!The xtilde values etc. are the real values of x, turned into dimensionless parameters via the 'thickness' d
! i.e. xtilde = x/d etc.

ztildei=-20
xtildei=-20
ztildef=20
xtildef=20
ztildestepfrac=0.1
xtildestepfrac=0.1

!The distance from the source to the interface parametised by d (in zspace)
dsourcetilde = 1

!Parameters: use command line args in future http://web.utah.edu/thorne/computing/Handy_Fortran_Tricks.pdf

eps1=1.0
mu1=1.0
eps2=10.0
mu2=1.0
ti = '45'


!angle of incidence (i.e. rotation of frame, but we are staying in non-rotated space)


PI=4.D0*DATAN(1.D0) 
!ensures maximum precision on any architechture apparently
thetai= (45.0/180.0)*PI 




n1=SQRT(eps1*mu1)
n2=SQRT(eps2*mu2)

if ((RealPart(n2) < 0 .and. RealPart(eps2) > 0) .or. &
 (RealPart(eps2) < 0 .and. mu2<0 .and. RealPart(n2) > 0 )) then
	n2 = -1 * n2
end if

xtildesize = anint(((xtildef-xtildei)/xtildestepfrac))
ztildesize = anint(((ztildef-ztildei)/ztildestepfrac))




!allocate arrays with desired shape
allocate(xtildearray(0:xtildesize))
allocate(ztildearray(0:ztildesize))
allocate(Eyrealarray(0:xtildesize,0:ztildesize))
allocate(Eykspacearray(0:xtildesize,0:ztildesize))


i = (0.0,1.0)
!Can just use normal functions as modern fortran can determine the type required, 
!specialist csqrt etc. are obsolete

eta = PI
!go from 0.1 to 5, dimensionless parameter equal to omega*d/c
! where d is the thickness of the slab and lambda is the free space wavelength of the incident light
! w and d and lambda are all replaced by eta


sigmatilde = 4.0
!sigma is form of gaussian beamwidth parameter, here taken to be the squared value, dimensions of area
!sigma tilde is sigma/d^2 so it is dimensionless parameter


tilen=LEN(TRIM(ti)) !this is just for filename purposes
kxtildeprimestepfrac = ((eta)/100.0) 


do p=0, xtildesize
	xtildearray(p)= xtildei + p*xtildestepfrac
end do
!print*, "x=", xtildesize," ", xtildearray(0)
do p=0, ztildesize
	ztildearray(p)= ztildei + p*ztildestepfrac
end do
!print*, "z=" ,ztildesize," ", ztildearray(0)
do p=0, 200
	kxtildeprimearray(p)= -(eta) + p*kxtildeprimestepfrac !used to be 0+ 
end do


do m=0, ztildesize
	do n=0, xtildesize
		Eykspacearray(n,m)=0
		
			do p=0, 200
			!avoid singularities
				kztildeprime=SQRT(eta**2 - kxtildeprimearray(p)**2)



				kxtilde=kxtildeprimearray(p)*cos(thetai) + kztildeprime*sin(thetai)

				kz1tilde=sqrt(eta**2 - kxtilde**2)
				kz2tilde=sqrt((n2*eta)**2 - kxtilde**2)



				if (aimag(kz2tilde) < 0.0005) then
					kz2tilde = real(kz2tilde)
				end if

				if (aimag(kz1tilde) < 0.0005) then
					kz1tilde = real(kz1tilde)
				end if

				if ((RealPart(kz2tilde) > 0) .and. (RealPart(n2) < 0)) then
					kz2tilde = -kz2tilde
				end if

				if (RealPart(kz1tilde) < 0) then
					kz1tilde = -kz1tilde
				end if


				chi = (kz2tilde*mu1)/(kz1tilde*mu2)

				A=exp(i*kz1tilde*dsourcetilde)
				B=exp(i*kz2tilde*dsourcetilde)
				C2=exp(i*kz2tilde*11*dsourcetilde)
				D=exp(i*kz1tilde*11*dsourcetilde)
				rnumerator=  ( 2*( ((chi + 1)* (C2**(-2))) + ((chi - 1)*(B**(-2))) ) )
				rdenominator= ( (((chi + 1)**2)*(C2**(-2))) - (((chi - 1)**2)*(B**(-2))) ) 
				r=  (rnumerator/rdenominator)  - 1
				if (chi==1) then
					!r=0
				end if 

				Ce= (2*A * (B**(-1)) * (C2**(-2)) * (chi + 1) ) / ( (((chi + 1)**2)*(C2**(-2)))  - (((chi - 1)**2)*(B**(-2))) ) 	
				De=(2*(chi - 1)*A*(B**(-1))) / ( (((chi + 1)**2)*(C2**(-2)))  - (((chi - 1)**2)*(B**(-2))) )
				t=(chi/D)*( (4*A* (B**(-1)) *(C2**(-1)))  / ( (((chi + 1)**2)*(C2**(-2)))  - (((chi - 1)**2)*(B**(-2))) ) )		



				etatest=eta**2-kxtilde**2

				if(ztildearray(m) <= dsourcetilde) then
					integral= ((1.0/sqrt(2*PI))*sqrt(sigmatilde)) &
					*(EXP( (-sigmatilde*((kxtilde*cos(thetai) - kz1tilde*sin(thetai) )**2)/2.0) & 
					 + (i* kxtilde * xtildearray(n)  ) + (i*kz1tilde*ztildearray(m) ) ) & 
			 		 *kxtildeprimestepfrac)

					integral2= ((1.0/sqrt(2*PI))*sqrt(sigmatilde)) &
					*(EXP( (-sigmatilde*((kxtilde*cos(thetai) - kz1tilde*sin(thetai) )**2)/2.0) & 
					 + (i* kxtilde * xtildearray(n)  ) + (-i*kz1tilde*ztildearray(m) ) ) &
			 		 *kxtildeprimestepfrac)


					Eykspacearray(n,m) = Eykspacearray(n,m) + (integral + r*integral2)
					if(m==0 .and. n==0) then
						!print*, "conserved= ",cdabs(r)**2+(kz2tilde/kz1tilde)*(cdabs(t)**2)
					end if

				elseif ( ztildearray(m) <= (11*dsourcetilde) ) then

					integral= ((1.0/sqrt(2*PI))*sqrt(sigmatilde)) &
					*(EXP( (-sigmatilde*((kxtilde*cos(thetai) - kz1tilde*sin(thetai) )**2)/2.0) & 
					 + (i* kxtilde * xtildearray(n)  ) + (i*kz2tilde*ztildearray(m) ) ) & 
			 		 *kxtildeprimestepfrac)

					integral2= ((1.0/sqrt(2*PI))*sqrt(sigmatilde)) &
					*(EXP( (-sigmatilde*((kxtilde*cos(thetai) - kz1tilde*sin(thetai) )**2)/2.0) & 
					 + (i* kxtilde * xtildearray(n)  ) + (-i*kz2tilde*ztildearray(m) ) ) &
			 		 *kxtildeprimestepfrac)


					Eykspacearray(n,m) = Eykspacearray(n,m) + (Ce*integral + De*integral2)
	


				else

					integral= ((1.0/sqrt(2*PI))*sqrt(sigmatilde)) &
					*(EXP( (-sigmatilde*((kxtilde*cos(thetai) - kz1tilde*sin(thetai) )**2)/2.0) & 
					 + (i* kxtilde* xtildearray(n)  ) + (i*kz1tilde*ztildearray(m) ) ) & 
			 		 *kxtildeprimestepfrac)


					Eykspacearray(n,m) = Eykspacearray(n,m) + (t*integral)
					if(m==0 .and. n==0) then
						!print*, " conserved= ",( cdabs(r)**2+(kz2tilde/kz1tilde)*(cdabs(t)**2) )
					end if
					



				end if
				
					if(etatest.LT.0.0d0) then
						Eykspacearray(n,m)=0
					end if

			end do 
	end do 
end do

write(filename,20) 'data/doubint',ti(1:tilen),'degs',eta,'eta', sigmatilde,'sigmatilde.dat'
open(unit=2,file= filename)
	


do m=0, ztildesize
	do n=0, xtildesize
		write(2,10) xtildearray(n), ztildearray(m), abs(Eykspacearray(n,m)), realpart(Eykspacearray(n,m))
	end do
end do
10	format(4e15.5,4e15.5,4e15.5,4e15.5)
20 	format(A,A,A,f3.1,A,f3.1,A)


end

!