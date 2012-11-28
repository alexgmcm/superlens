!resolutionsuperlenscli.f90
!Fixing major bug due to incorrect integral limits"
!therefore change kxtilde for kxtildeprime in integral
!Single interface code using rotated frame approach but translated
!to normal frame.
!Supports command line interface for variables and calls to MATLAB for plotting.
!add evanescent components
!COPY FOR WORKING ON

!have CLI just for second interface (to vary slab thickness) and thetamax, the maximum incident angle assumed to hit the slab
!CLI: ./test.out SECONDINTERFACE THETAMAX

!
implicit none !doesn't assume types from names, must be declared explicitly
!use strings
! Declare stuff here - check these are all necessary
double precision :: thetai, eta, xtildestepfrac, ztildestepfrac, kxtildeprimestepfrac, dsourcetilde, secondinterface
double precision :: xtildei, ztildei, c=3D8, xtildef, ztildef, PI, eps1, mu1, mu2, etatest, sigmatilde, smallval, thetamax 
double precision :: thetamaxrad, cutkxtildeprimestepfrac
complex*16 :: A, B, C2, D, i, integral, integral2, chi, n1, n2, eps2
complex*16 :: r, t, kztildeprime, kz1tilde,kz2tilde,kxtilde, Ce, De, rnumerator, rdenominator
complex*16 :: cutkztildeprime, cutkxtilde, cutkz1tilde, cutkz2tilde, cutchi, cutA, cutB, cutC2, cutD, cutrnumerator, cutrdenominator
complex*16 :: cutCe, cutDe, cutt, cutr
integer*4 :: m, n, p, xtildesize, ztildesize,tilen, errflag
character :: filename*150, ti*10, cmd1*50, cmd2*50 
double precision, dimension(:), allocatable :: xtildearray
double precision, dimension(:), allocatable :: ztildearray
complex*16, dimension(0:200) :: kxtildeprimearray
complex*16, dimension(0:200) :: cutkxtildeprimearray
complex*16, dimension(:,:), allocatable :: Eyrealarray
complex*16, dimension(:,:), allocatable :: Eykspacearray





!The xtilde values etc. are the real values of x, turned into dimensionless parameters via the 'thickness' d
! i.e. xtilde = x/d etc.

ztildei=-5
xtildei=-10
ztildef=15
xtildef=10
ztildestepfrac=0.1
xtildestepfrac=0.01

!The distance from the source to the interface parametised by d (in zspace)
dsourcetilde = 1

!All distances are essentially in terms of dsourcetilde, however they will simply be given as numbers
!as this makes passing them as commmand line arguments far far easier (and is equivalent as dsourcetilde=1)

!position of second interface
!secondinterface = 3*dsourcetilde

!Parameters: use command line args in future http://web.utah.edu/thorne/computing/Handy_Fortran_Tricks.pdf

eps1=1.0
mu1=1.0
eps2=-1.0
mu2=-1.0
ti = '0' !also change thetai

!Remember negative refraction is when BOTH eps2 AND mu2 are negative, not just eps2!!!
!change thetai AND ti
!angle of incidence (i.e. rotation of frame, but we are staying in non-rotated space)


PI=4.D0*DATAN(1.D0) 
!ensures maximum precision on any architechture apparently
smallval = 0.1
thetai= (smallval/180.0)*PI !also change ti

CALL GETARG(2,cmd2)
READ(UNIT=cmd2, FMT=*) thetamax
print*, "thetamax=", thetamax
thetamaxrad = (thetamax/180.0)*PI



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

sigmatilde=0.001
CALL GETARG(1,cmd1)
READ(UNIT=cmd1, FMT=*) secondinterface
!CALL VALUE(cmd, sigmatilde, errflag) 
!print*, "errflag=", errflag
!sigma is form of gaussian beamwidth parameter, here taken to be the squared value, dimensions of area
!sigma tilde is sigma/d^2 so it is dimensionless parameter
print*, "secondinterface=", secondinterface

tilen=LEN(TRIM(ti)) !this is just for filename purposes
cutkxtildeprimestepfrac = ((1*eta*SIN(thetamaxrad))/100.0) 
kxtildeprimestepfrac = ((1*eta)/100.0)

do p=0, xtildesize
	xtildearray(p)= xtildei + p*xtildestepfrac
end do
!print*, "x=", xtildesize," ", xtildearray(0)
do p=0, ztildesize
	ztildearray(p)= ztildei + p*ztildestepfrac
end do
!print*, "z=" ,ztildesize," ", ztildearray(0)
do p=0, 200
	kxtildeprimearray(p)= -(1*eta) + p*kxtildeprimestepfrac
	cutkxtildeprimearray(p)= -(1*eta*SIN(thetamaxrad)) + p*cutkxtildeprimestepfrac !>1*eta corresponds to including dark modes 
end do


do m=0, ztildesize
	do n=0, xtildesize
		Eykspacearray(n,m)=0
		
			do p=0, 200
			!avoid singularities
				kztildeprime=SQRT(eta**2 - kxtildeprimearray(p)**2)
				cutkztildeprime=SQRT(eta**2 - cutkxtildeprimearray(p)**2)


				kxtilde=kxtildeprimearray(p)*cos(thetai) + kztildeprime*sin(thetai)
				cutkxtilde=cutkxtildeprimearray(p)*cos(thetai) + cutkztildeprime*sin(thetai)

				kz1tilde=sqrt(eta**2 - kxtilde**2)
				cutkz1tilde=sqrt(eta**2 - cutkxtilde**2)

				kz2tilde=sqrt((n2*eta)**2 - kxtilde**2)
				cutkz2tilde=sqrt((n2*eta)**2 - cutkxtilde**2)

				if (aimag(kz2tilde)<0) then
					kz2tilde = cmplx(-real(kz2tilde), -aimag(kz2tilde))

				end if

				if (aimag(kz1tilde)<0) then
					kz1tilde = cmplx(-real(kz1tilde), -aimag(kz1tilde))
				end if

				if (aimag(kxtilde)<0) then
					kxtilde = cmplx(-real(kxtilde), -aimag(kxtilde))
				end if


				if (aimag(cutkz2tilde)<0) then
					cutkz2tilde = cmplx(-real(cutkz2tilde), -aimag(cutkz2tilde))

				end if

				if (aimag(cutkz1tilde)<0) then
					cutkz1tilde = cmplx(-real(cutkz1tilde), -aimag(cutkz1tilde))
				end if

				if (aimag(cutkxtilde)<0) then
					cutkxtilde = cmplx(-real(cutkxtilde), -aimag(cutkxtilde))
				end if



				chi = (kz2tilde*mu1)/(kz1tilde*mu2)
				A=exp(i*kz1tilde*dsourcetilde)
				B=exp(i*kz2tilde*dsourcetilde)
				C2=exp(i*kz2tilde*secondinterface)
				D=exp(i*kz1tilde*secondinterface)
				rnumerator=  ( 2*( ((chi + 1)* (C2**(-2))) + ((chi - 1)*(B**(-2))) ) )
				rdenominator= ( (((chi + 1)**2)*(C2**(-2))) - (((chi - 1)**2)*(B**(-2))) ) 
				r=  (rnumerator/rdenominator)  - 1
				


				Ce= (2*A * (B**(-1)) * (C2**(-2)) * (chi + 1) ) / ( (((chi + 1)**2)*(C2**(-2)))  - (((chi - 1)**2)*(B**(-2))) ) 	
				De=(2*(chi - 1)*A*(B**(-1))) / ( (((chi + 1)**2)*(C2**(-2)))  - (((chi - 1)**2)*(B**(-2))) )
				t=(chi/D)*( (4*A* (B**(-1)) *(C2**(-1)))  / ( (((chi + 1)**2)*(C2**(-2)))  - (((chi - 1)**2)*(B**(-2))) ) )		



				cutchi = (cutkz2tilde*mu1)/(cutkz1tilde*mu2)
				cutA=exp(i*cutkz1tilde*dsourcetilde)
				cutB=exp(i*cutkz2tilde*dsourcetilde)
				cutC2=exp(i*cutkz2tilde*secondinterface)
				cutD=exp(i*cutkz1tilde*secondinterface)
				cutrnumerator=  ( 2*( ((cutchi + 1)* (cutC2**(-2))) + ((cutchi - 1)*(cutB**(-2))) ) )
				cutrdenominator= ( (((cutchi + 1)**2)*(cutC2**(-2))) - (((cutchi - 1)**2)*(cutB**(-2))) ) 
				cutr=  (cutrnumerator/cutrdenominator)  - 1
				


				cutCe= (2*cutA * (cutB**(-1)) * (cutC2**(-2)) * (cutchi + 1) ) & 
					/ ( (((cutchi + 1)**2)*(cutC2**(-2)))  - (((cutchi - 1)**2)*(cutB**(-2))) ) 	
				cutDe=(2*(cutchi - 1)*cutA*(cutB**(-1))) / ( (((cutchi + 1)**2)*(cutC2**(-2))) - (((cutchi - 1)**2)*(cutB**(-2))))
				cutt=(cutchi/cutD)*( (4*cutA* (cutB**(-1)) *(cutC2**(-1))) &
					 / ( (((cutchi + 1)**2)*(cutC2**(-2)))  - (((cutchi - 1)**2)*(cutB**(-2))) ) )	



				etatest=eta**2-kxtilde**2

				if(ztildearray(m) <= dsourcetilde) then
					integral= ((1.0/sqrt(2*PI))*sqrt(sigmatilde)) &
					*(EXP( (-sigmatilde*(( kxtildeprimearray(p) )**2)/2.0) & 
					 + (i* kxtilde * xtildearray(n)  ) + (i*kz1tilde*ztildearray(m) ) ) & 
			 		 *kxtildeprimestepfrac)

					integral2= ((1.0/sqrt(2*PI))*sqrt(sigmatilde)) &
					*(EXP( (-sigmatilde*(( cutkxtildeprimearray(p) )**2)/2.0) & 
					 + (i* cutkxtilde * xtildearray(n)  ) + (-i*cutkz1tilde*ztildearray(m) ) ) &
			 		 *cutkxtildeprimestepfrac)


					Eykspacearray(n,m) = Eykspacearray(n,m) + (integral + cutr*integral2)
					if(m==0 .and. n==0) then
						!print*, "conserved= ",cdabs(r)**2+(kz2tilde/kz1tilde)*(cdabs(t)**2)
					end if

				elseif ( ztildearray(m) <= (secondinterface) ) then

					integral= ((1.0/sqrt(2*PI))*sqrt(sigmatilde)) &
					*(EXP( (-sigmatilde*(( cutkxtildeprimearray(p) )**2)/2.0) & 
					 + (i* cutkxtilde * xtildearray(n)  ) + (i*cutkz2tilde*ztildearray(m) ) ) & 
			 		 *cutkxtildeprimestepfrac)

					integral2= ((1.0/sqrt(2*PI))*sqrt(sigmatilde)) &
					*(EXP( (-sigmatilde*(( cutkxtildeprimearray(p) )**2)/2.0) & 
					 + (i* cutkxtilde * xtildearray(n)  ) + (-i*cutkz2tilde*ztildearray(m) ) ) &
			 		 *cutkxtildeprimestepfrac)


					Eykspacearray(n,m) = Eykspacearray(n,m) + (cutCe*integral + cutDe*integral2)
	


				else

					integral= ((1.0/sqrt(2*PI))*sqrt(sigmatilde)) &
					*(EXP( (-sigmatilde*(( cutkxtildeprimearray(p) )**2)/2.0) & 
					 + (i* cutkxtilde* xtildearray(n)  ) + (i*cutkz1tilde*ztildearray(m) ) ) & 
			 		 *cutkxtildeprimestepfrac)


					Eykspacearray(n,m) = Eykspacearray(n,m) + (cutt*integral)
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

write(filename,20) 'data/res',thetamax,'degs',eta,'eta', sigmatilde,'sigmatilde',secondinterface,'secint.dat'
open(unit=2,file= filename)
	


do m=0, ztildesize
	do n=0, xtildesize
		write(2,10) xtildearray(n), ztildearray(m), (abs(Eykspacearray(n,m))**2), realpart(Eykspacearray(n,m)) !intensity so square
	end do
end do
10	format(4e15.5,4e15.5,4e15.5,4e15.5)
20 	format(A,f4.1,A,f3.1,A,f5.3,A,f3.1,A)

!cmd='./matlab_batcher.sh superlenscliscript ', sigmatilde
!write (cmd, "(A39,I2)") "hello", 10

!CALL SYSTEM(cmd)

end

!