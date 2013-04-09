!Supports command line interface for variables and calls to MATLAB for plotting.
!added evanescent components
!Command Line Interface is used to set various values and determine which code is run
!The code
!CLI EXAMPLE: ./test.out SECONDINTERFACE THETAMAX ETALIMIT MODEFLAG IMAGEFLAG SIGMATILDE
PROGRAM masterresolutionsuperlenscli
implicit none !doesn't assume types from names, must be declared explicitly
! Declare variables here:
double precision :: thetai, eta, xtildestepfrac, ztildestepfrac, kxtildestepfrac, dsourcetilde, secondinterface
double precision :: xtildei, ztildei, c=3D8, xtildef, ztildef, PI, eps1, mu1, mu2, etatest, sigmatilde, smallval, thetamax 
double precision :: thetamaxrad, ztilde, xtilde
complex*16 :: A, B, C2, D, i, integral, integral2, chi, n1, n2, eps2, addterm
complex*16 :: r, t, kz1tilde,kz2tilde,kxtilde, Ce, De, rnumerator, rdenominator, Eykspace
integer*4 :: m, n, p, xtildesize, ztildesize,tilen, errflag, numkxpoints,darknumkxpoints,cutnumkxpoints, etalimit, modeflag
integer*4 :: imageflag
character :: filename*150, ti*10, cmd1*50, cmd2*50, cmd3*2, etalimitstring*2, cmd4*1, cmd5*1, imagestring*20, modestring*20
character :: kxfilename*150, cmd6*10
complex*16 :: darkkz1tilde, darkkz2tilde, darkdenominator, darkA, darkC2, darkD, darkT




PI=4.D0*DATAN(1.D0) 
!ensures maximum precision on any architechture


CALL GETARG(1,cmd1)
READ(UNIT=cmd1, FMT=*) secondinterface
!position of second interface

CALL GETARG(2,cmd2)
READ(UNIT=cmd2, FMT=*) thetamax
print*, "thetamax=", thetamax
thetamaxrad = (thetamax/180.0)*PI
!thetacutoff value entered in degrees so convert to radians


CALL GETARG(3,cmd3)
READ(UNIT=cmd3, FMT=*) etalimit

CALL GETARG(4,cmd4)
READ(UNIT=cmd4, FMT=*) modeflag
!modeflag=0 proponly, 1 combined, 2 darkonly


CALL GETARG(5,cmd5)
READ(UNIT=cmd5, FMT=*) imageflag
!imageflag=0 full 2d, 1 image, 2 source

CALL GETARG(6,cmd6)
READ(UNIT=cmd6, FMT=*) sigmatilde
!sigma is form of gaussian beamwidth parameter, here taken to be the squared value, dimensions of area
!sigma tilde is sigma/d^2 so it is dimensionless parameter

!The xtilde values etc. are the real values of x, turned into dimensionless parameters via 
!the source to first interface distance d_s
! i.e. xtilde = x/d_s etc.
xtildei=-1
xtildef=1

if(imageflag==0) then
	xtildei=-10
	xtildef=10
end if


ztildei=0
ztildef=10

ztildestepfrac=0.1
xtildestepfrac=0.001
if(imageflag==0) then
	xtildestepfrac=0.01
end if

!The distance from the source to the interface parametised by d (in zspace)
dsourcetilde = 1

!All distances are essentially in terms of dsourcetilde, however they will simply be given as numbers
!as this makes passing them as commmand line arguments far far easier (and is equivalent as dsourcetilde=1)

eps1=1.0
mu1=1.0
eps2=-1.0
mu2=-1.0
!Remember negative refraction is when BOTH eps2 AND mu2 are negative, not just eps2!!!

n1=SQRT(eps1*mu1)
n2=SQRT(eps2*mu2)

if ((RealPart(n2) < 0 .and. RealPart(eps2) > 0) .or. &
	(RealPart(eps2) < 0 .and. mu2<0 .and. RealPart(n2) > 0 )) then
n2 = -1 * n2
end if

xtildesize = anint(((xtildef-xtildei)/xtildestepfrac))
ztildesize = anint(((ztildef-ztildei)/ztildestepfrac))


i = (0.0,1.0)
!Can just use normal functions as modern fortran can determine the type required, 
!specialist csqrt etc. are obsolete

eta = PI
!dimensionless parameter equal to omega*d_s/c

print*, "secondinterface=", secondinterface

tilen=LEN(TRIM(ti)) !this is just for filename purposes

numkxpoints = 200
kxtildestepfrac = ((2.0*eta)/numkxpoints)
darknumkxpoints = anint(((etalimit-1)*eta)/kxtildestepfrac)
cutnumkxpoints = anint((2.0*eta*SIN(thetamaxrad))/kxtildestepfrac)

SELECT CASE (modeflag)
	CASE (0)
		modestring='proponly'
	CASE (1)
		modestring='combined'
	CASE (2)
		modestring='darkonly'
	CASE DEFAULT
		print*, "modeflag ERROR: ", modeflag
END SELECT 


SELECT CASE	(imageflag)
	CASE (0)
		imagestring='2D'
	CASE (1)
		imagestring='imageline'
	CASE (2)
		imagestring='sourceline'
	CASE DEFAULT
		print*, "imageflag ERROR: ", imageflag
END SELECT

write(etalimitstring, 30) etalimit
write(filename,20) 'data/',trim(adjustl(modestring)),trim(adjustl(imagestring)),thetamax,'degs',eta,'eta', &
sigmatilde,'sigmatilde',secondinterface, 'secint',trim(adjustl(etalimitstring)),'etalimit.dat'

open(unit=2,file= filename)

do m=0, ztildesize
	ztilde= ztildei + (m*ztildestepfrac)

	if(imageflag==1) then
		ztilde=4.0
	end if

	if(imageflag==2) then
		ztilde = 0.0
	end if




	do n=0, xtildesize
		xtilde= xtildei + (n*xtildestepfrac)
		Eykspace=0
		!REGION ONE:
		if(ztilde <= dsourcetilde) then

			if(modeflag/=0) then
			!DARK PARTS, between -infty and -eta, and eta and infty, sum the two, infty is taken as 5*eta
			do p=1, darknumkxpoints
				kxtilde= -(etalimit*eta) + p*kxtildestepfrac !this part depends on limits
				call SHAREDINTEGRALCODE()
				integral=(exp(-darkkz1tilde*ztilde)+ &
					darkA*exp(darkkz1tilde*(ztilde-dsourcetilde)) ) * kxtildestepfrac *exp(i*kxtilde*xtilde) *&
				((1.0/sqrt(2*PI))*sqrt(sigmatilde)) *(EXP( (-sigmatilde*(( kxtilde )**2)/2.0)))
				Eykspace = Eykspace + integral
				! integral code will depend on limits/area so can't go in subroutine
				
				
			end do

			!Dark parts upper limits
			do p=1, darknumkxpoints
				kxtilde= (eta) + p*kxtildestepfrac !this part depends on limits
				call SHAREDINTEGRALCODE()
				integral=(exp(-darkkz1tilde*ztilde)+ &
					darkA*exp(darkkz1tilde*(ztilde-dsourcetilde)) ) * kxtildestepfrac * exp(i*kxtilde*xtilde)*&
				((1.0/sqrt(2*PI))*sqrt(sigmatilde)) *(EXP( (-sigmatilde*(( kxtilde )**2)/2.0)))
				Eykspace = Eykspace + integral
			end do
		end if

		if(modeflag /= 2) then
			!light part non-truncated (incident wave)

			do p=1, numkxpoints
				kxtilde= (-eta) + p*kxtildestepfrac !this part depends on limits
				call SHAREDINTEGRALCODE()
				integral= ((1.0/sqrt(2*PI))*sqrt(sigmatilde)) &
				*(EXP( (-sigmatilde*(( kxtilde )**2)/2.0) & 
					+ (i* kxtilde * xtilde  ) + (i*kz1tilde*ztilde ) ) & 
				*kxtildestepfrac)
				Eykspace = Eykspace +integral
			end do
			!light part truncated (reflected wave)
			do p=1, cutnumkxpoints
				kxtilde= (-eta*SIN(thetamaxrad)) + p*kxtildestepfrac !this part depends on limits
				call SHAREDINTEGRALCODE()
				integral= ((1.0/sqrt(2*PI))*sqrt(sigmatilde)) &
				*(EXP( (-sigmatilde*(( kxtilde )**2)/2.0) & 
					+ (i* kxtilde * xtilde  ) + (-i*kz1tilde*ztilde ) ) & 
				*kxtildestepfrac)
				Eykspace = Eykspace + (r*integral)
			end do
		end if
			!REGION TWO
			elseif ( ztilde <= (secondinterface) ) then

				if(modeflag/=0) then
			!DARK PARTS, between -infty and -eta, and eta and infty, sum the two, infty is taken as 5*eta
			do p=1, darknumkxpoints
				kxtilde= -(etalimit*eta) + p*kxtildestepfrac !this part depends on limits
				call SHAREDINTEGRALCODE()
				integral=( darkC2*exp(-darkkz2tilde*(ztilde-dsourcetilde)) +&
					darkD*exp(-darkkz2tilde* (3*dsourcetilde - ztilde) ) ) * kxtildestepfrac *exp(i*kxtilde*xtilde)*&
				((1.0/sqrt(2*PI))*sqrt(sigmatilde)) *(EXP( (-sigmatilde*(( kxtilde )**2)/2.0)))
				Eykspace = Eykspace + integral
			end do

			!Dark parts upper limits
			do p=1, darknumkxpoints
				kxtilde= (eta) + p*kxtildestepfrac !this part depends on limits
				call SHAREDINTEGRALCODE()
				integral=( darkC2*exp(-darkkz2tilde*(ztilde-dsourcetilde)) +&
					darkD*exp(-darkkz2tilde* (3*dsourcetilde - ztilde) ) ) * kxtildestepfrac *exp(i*kxtilde*xtilde)*&
				((1.0/sqrt(2*PI))*sqrt(sigmatilde)) *(EXP( (-sigmatilde*(( kxtilde )**2)/2.0)))
				Eykspace = Eykspace + integral
			end do
		end if
		if(modeflag/=2) then
			!light part truncated (all waves)
			do p=1, cutnumkxpoints
				kxtilde= (-eta*SIN(thetamaxrad)) + p*kxtildestepfrac !this part depends on limits
				call SHAREDINTEGRALCODE()
				integral= ((1.0/sqrt(2*PI))*sqrt(sigmatilde)) &
				*(EXP( (-sigmatilde*(( kxtilde )**2)/2.0) & 
					+ (i* kxtilde * xtilde  ) + (i*kz2tilde*ztilde ) ) & 
				*kxtildestepfrac)

				integral2= ((1.0/sqrt(2*PI))*sqrt(sigmatilde)) &
				*(EXP( (-sigmatilde*(( kxtilde )**2)/2.0) & 
					+ (i* kxtilde * xtilde  ) + (-i*kz2tilde*ztilde ) ) &
				*kxtildestepfrac)
				Eykspace = Eykspace + (Ce*integral + De*integral2)
			end do
		end if
			!REGION THREE
			else
				if(modeflag/=0) then
			!DARK PARTS, between -infty and -eta, and eta and infty, sum the two, infty is taken as 5*eta
			do p=1, darknumkxpoints
				kxtilde= -(etalimit*eta) + p*kxtildestepfrac !this part depends on limits
				call SHAREDINTEGRALCODE()
				integral=(darkT*exp(-darkkz1tilde*(ztilde- 3*dsourcetilde))) * kxtildestepfrac *exp(i*kxtilde*xtilde)*&
				((1.0/sqrt(2*PI))*sqrt(sigmatilde)) *(EXP( (-sigmatilde*(( kxtilde )**2)/2.0)))
				Eykspace = Eykspace + integral

			end do
			!Dark parts upper limits
			do p=1, darknumkxpoints
				kxtilde= (eta) + p*kxtildestepfrac !this part depends on limits
				call SHAREDINTEGRALCODE()
				integral=(darkT*exp(-darkkz1tilde*(ztilde- 3*dsourcetilde))) * kxtildestepfrac *exp(i*kxtilde*xtilde)*&
				((1.0/sqrt(2*PI))*sqrt(sigmatilde)) *(EXP( (-sigmatilde*(( kxtilde )**2)/2.0)))
				Eykspace = Eykspace + integral

			end do
		end if
		if(modeflag/=2) then
			!light part truncated (transmitted wave)
			do p=1, cutnumkxpoints
				kxtilde= (-eta*SIN(thetamaxrad)) + p*kxtildestepfrac !this part depends on limits
				call SHAREDINTEGRALCODE()
				integral= ((1.0/sqrt(2*PI))*sqrt(sigmatilde))
				integral=integral*(EXP( (-sigmatilde*(( kxtilde )**2)/2.0)))
				integral=integral*(EXP(i*kxtilde*xtilde))
				integral=integral*(EXP((i* kxtilde* xtilde )))
				integral=integral*(EXP(i*kz1tilde*ztilde))
				integral = integral*kxtildestepfrac
				addterm = (t*integral)
				Eykspace = Eykspace + addterm
			end do
		end if
	end if
	
		write(2,40) xtilde, ztilde, (abs(Eykspace))**2, abs(Eykspace) !intensity - therefore squared
	end do

	if (imageflag /= 0) then
		!if we only need one line of data then quit after getting it
		EXIT
	end if

end do


40	format(4e15.5,4e15.5,4e15.5,4e15.5)
15  format(4e15.5,4e15.5,4e15.5,4e15.5,4e15.5)
20 	format(A,A,A,f4.1,A,f3.1,A,f5.3,A,f3.1,A,A,A)
30  format(I2)

CONTAINS
SUBROUTINE SHAREDINTEGRALCODE()

kz1tilde=sqrt(eta**2 - kxtilde**2)
kz2tilde=sqrt((n2*eta)**2 - kxtilde**2)
		!TO MAKE SURE THAT  i ISNT TAKEN INTO THE KZ PARTS FOR DARK WAVES!!!!
		darkkz1tilde=sqrt(abs((n1*eta)**2 - kxtilde**2))
		darkkz2tilde=sqrt(abs((n2*eta)**2 - kxtilde**2))

		if (aimag(kz2tilde)<0) then
			kz2tilde = cmplx(-real(kz2tilde), -aimag(kz2tilde))
		end if
		if (aimag(kz1tilde)<0) then
			kz1tilde = cmplx(-real(kz1tilde), -aimag(kz1tilde))
		end if
		if (aimag(kxtilde)<0) then
			kxtilde = cmplx(-real(kxtilde), -aimag(kxtilde))
		end if
		!makes sure that the complex square root is handled correctly
		!no exponentially growing parts
		chi = (kz2tilde*mu1)/(kz1tilde*mu2)

		A=exp(i*kz1tilde*dsourcetilde)
		B=exp(i*kz2tilde*dsourcetilde)
		C2=exp(i*kz2tilde*secondinterface)
		D=exp(i*kz1tilde*secondinterface)
		rnumerator=  ( 2*( ((chi + 1)* (C2**(-2))) + ((chi - 1)*(B**(-2))) ) )
		rdenominator= ( (((chi + 1)**2)*(C2**(-2))) - (((chi - 1)**2)*(B**(-2))) ) 
		r= (rnumerator/rdenominator)  - 1

		Ce= (2*A * (B**(-1)) * (C2**(-2)) * (chi + 1) ) / ( (((chi + 1)**2)*(C2**(-2)))  - (((chi - 1)**2)*(B**(-2))) ) 	
		De=(2*(chi - 1)*A*(B**(-1))) / ( (((chi + 1)**2)*(C2**(-2)))  - (((chi - 1)**2)*(B**(-2))) )
		t=(chi/D)*( (4*A* (B**(-1)) *(C2**(-1)))  / ( (((chi + 1)**2)*(C2**(-2)))  - (((chi - 1)**2)*(B**(-2))) ) )		

		darkdenominator=(((chi+1)**2)-((chi-1)**2)*exp(-4*darkkz2tilde*dsourcetilde)) 
		darkD=((2*(chi-1)*exp(-darkkz1tilde*dsourcetilde)*exp(-2*darkkz2tilde)*dsourcetilde))/darkdenominator
		darkT=(4*chi*exp(-darkkz1tilde*dsourcetilde)*exp(-2*darkkz2tilde*dsourcetilde))/darkdenominator
		darkC2=(2*exp(-darkkz1tilde*dsourcetilde)*(chi+1))/darkdenominator
		!darkA=darkC2 + darkD*exp(-2*darkkz2tilde*dsourcetilde) - exp(-darkkz1tilde*dsourcetilde)
		darkA=0 !should be zero but floating point is causing problems, so set as 0

		RETURN
	END SUBROUTINE

END PROGRAM