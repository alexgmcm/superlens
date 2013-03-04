!darkresolutionsuperlenscli.f90
!Fixing major bug due to incorrect integral limits"
!therefore change kxtilde for kxtilde in integral
!Single interface code using rotated frame approach but translated
!to normal frame.
!Supports command line interface for variables and calls to MATLAB for plotting.
!added evanescent components
!COPY FOR WORKING ON

!have CLI just for second interface (to vary slab thickness) and thetamax, the maximum incident angle assumed to hit the slab
!CLI: ./test.out SECONDINTERFACE THETAMAX
PROGRAM masterresolutionsuperlenscli
implicit none !doesn't assume types from names, must be declared explicitly
!use strings
! Declare stuff here - check these are all necessary
double precision :: thetai, eta, xtildestepfrac, ztildestepfrac, kxtildestepfrac, dsourcetilde, secondinterface
double precision :: xtildei, ztildei, c=3D8, xtildef, ztildef, PI, eps1, mu1, mu2, etatest, sigmatilde, smallval, thetamax 
double precision :: thetamaxrad, ztilde, xtilde
complex*16 :: A, B, C2, D, i, integral, integral2, chi, n1, n2, eps2, addterm
complex*16 :: r, t, kz1tilde,kz2tilde,kxtilde, Ce, De, rnumerator, rdenominator, Eykspace
integer*4 :: m, n, p, xtildesize, ztildesize,tilen, errflag, numkxpoints,darknumkxpoints,cutnumkxpoints, etalimit, modeflag
integer*4 :: imageflag
character :: filename*150, ti*10, cmd1*50, cmd2*50, cmd3*2, etalimitstring*2, cmd4*1, cmd5*1, imagestring*20, modestring*20
character :: kxfilename*150
complex*16 :: darkkz1tilde, darkkz2tilde, darkdenominator, darkA, darkC2, darkD, darkT




PI=4.D0*DATAN(1.D0) 
!ensures maximum precision on any architechture apparently


CALL GETARG(1,cmd1)
READ(UNIT=cmd1, FMT=*) secondinterface

CALL GETARG(2,cmd2)
READ(UNIT=cmd2, FMT=*) thetamax
print*, "thetamax=", thetamax
thetamaxrad = (thetamax/180.0)*PI

CALL GETARG(3,cmd3)
READ(UNIT=cmd3, FMT=*) etalimit

CALL GETARG(4,cmd4)
READ(UNIT=cmd4, FMT=*) modeflag
!modeflag=0 proponly, 1 combined, 2 darkonly


CALL GETARG(5,cmd5)
READ(UNIT=cmd5, FMT=*) imageflag
!imageflag=0 full 2d, 1 image, 2 source

!The xtilde values etc. are the real values of x, turned into dimensionless parameters via the 'thickness' d
! i.e. xtilde = x/d etc.
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

!position of second interface
!secondinterface = 3*dsourcetilde

!Parameters: use command line args in future http://web.utah.edu/thorne/computing/Handy_Fortran_Tricks.pdf

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
!go from 0.1 to 5, dimensionless parameter equal to omega*d/c
! where d is the thickness of the slab and lambda is the free space wavelength of the incident light
! w and d and lambda are all replaced by eta

sigmatilde=0.001

!CALL VALUE(cmd, sigmatilde, errflag) 
!print*, "errflag=", errflag
!sigma is form of gaussian beamwidth parameter, here taken to be the squared value, dimensions of area
!sigma tilde is sigma/d^2 so it is dimensionless parameter
print*, "secondinterface=", secondinterface

tilen=LEN(TRIM(ti)) !this is just for filename purposes

numkxpoints = 200
!cutkxtildestepfrac = ((2*eta*SIN(thetamaxrad))/numkxpoints) 
kxtildestepfrac = ((2.0*eta)/numkxpoints)
darknumkxpoints = anint(((etalimit-1)*eta)/kxtildestepfrac)
cutnumkxpoints = anint((2.0*eta*SIN(thetamaxrad))/kxtildestepfrac)




!etalimit=3

!darkstepfrac= (((etalimit-1)*eta)/numkxpoints)
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

write(kxfilename,20) 'data/kx',trim(adjustl(modestring)),trim(adjustl(imagestring)),thetamax,'degs',eta,'eta', &
sigmatilde,'sigmatilde',secondinterface, 'secint',trim(adjustl(etalimitstring)),'etalimit.dat'

open(unit=2,file= filename)
open(unit=3,file= kxfilename)
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
				!etatest=eta**2-kxtilde**2 ! ????? what does this do?
				! integral code will depend on limits/area so can't go in subroutine
				! need to check what form the equation is for dark components
				!code above here below subroutine part can be moved into subroutine though?
				!try using plane waves for dark parts, turn off light modes (comment out)
				
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
			!print*, "Doing dark parts. R3"
			do p=1, darknumkxpoints
				kxtilde= -(etalimit*eta) + p*kxtildestepfrac !this part depends on limits
				call SHAREDINTEGRALCODE()
				integral=(darkT*exp(-darkkz1tilde*(ztilde- 3*dsourcetilde))) * kxtildestepfrac *exp(i*kxtilde*xtilde)*&
				((1.0/sqrt(2*PI))*sqrt(sigmatilde)) *(EXP( (-sigmatilde*(( kxtilde )**2)/2.0)))
				Eykspace = Eykspace + integral

				if (ztilde==4.0 .and. xtilde<1.0E-05 .and. xtilde>-1.0E-05) then
					write(3,15) abs(kxtilde), abs(darkA), abs(darkC2), abs(darkD), abs(darkT)
					!imag part is always 0 anyway
				end if

			end do
			!Dark parts upper limits
			do p=1, darknumkxpoints
				kxtilde= (eta) + p*kxtildestepfrac !this part depends on limits
				call SHAREDINTEGRALCODE()
				integral=(darkT*exp(-darkkz1tilde*(ztilde- 3*dsourcetilde))) * kxtildestepfrac *exp(i*kxtilde*xtilde)*&
				((1.0/sqrt(2*PI))*sqrt(sigmatilde)) *(EXP( (-sigmatilde*(( kxtilde )**2)/2.0)))
				Eykspace = Eykspace + integral

				if (ztilde==4.0 .and. xtilde<1.0E-05 .and. xtilde>-1.0E-05) then
					write(3,15) abs(kxtilde), abs(darkA), abs(darkC2), abs(darkD), abs(darkT)
				end if

			end do
		end if
		if(modeflag/=2) then
			!light part truncated (transmitted wave)
			!print*, "Doing light parts. R3"
			do p=1, cutnumkxpoints
				kxtilde= (-eta*SIN(thetamaxrad)) + p*kxtildestepfrac !this part depends on limits
				call SHAREDINTEGRALCODE()
				integral= ((1.0/sqrt(2*PI))*sqrt(sigmatilde))
				!print*,"1 ", integral
				integral=integral*(EXP( (-sigmatilde*(( kxtilde )**2)/2.0)))
				!print*,"2 ", integral
				integral=integral*(EXP(i*kxtilde*xtilde))
				!print*,"3 ", integral
				integral=integral*(EXP((i* kxtilde* xtilde )))
				!print*,"4 ", integral
				integral=integral*(EXP(i*kz1tilde*ztilde))
				!print*,"5 ", integral
				integral = integral*kxtildestepfrac
				!print*,"6 ", integral
				!print*,"t ", t
				addterm = (t*integral)
				!print*,"addterm ", addterm
				!print*,"Eykspace: ", Eykspace
				Eykspace = Eykspace + addterm
				!print*,"Eykspace: ", Eykspace
				
				if (ztilde==4.0 .and. xtilde<1.0E-05 .and. xtilde>-1.0E-05) then
					write(3,15) abs(kxtilde), abs(r), abs(Ce), abs(De), abs(t)
				end if

			end do
		end if
	end if
	
		write(2,40) xtilde, ztilde, (abs(Eykspace))**2, abs(Eykspace) !intensity - therefore squared
	end do

	if (imageflag /= 0) then
		EXIT
	end if

end do


40	format(4e15.5,4e15.5,4e15.5,4e15.5)
15  format(4e15.5,4e15.5,4e15.5,4e15.5,4e15.5)
20 	format(A,A,A,f4.1,A,f3.1,A,f5.3,A,f3.1,A,A,A)
30  format(I2)
!40 	format()

!cmd='./matlab_batcher.sh superlenscliscript ', sigmatilde
!write (cmd, "(A39,I2)") "hello", 10

!CALL SYSTEM(cmd)



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

		!print*, "r=", r, " ce=",Ce," De=",De," t=",t

		darkdenominator=(((chi+1)**2)-((chi-1)**2)*exp(-4*darkkz2tilde*dsourcetilde)) 
		darkD=((2*(chi-1)*exp(-darkkz1tilde*dsourcetilde)*exp(-2*darkkz2tilde)*dsourcetilde))/darkdenominator
		darkT=(4*chi*exp(-darkkz1tilde*dsourcetilde)*exp(-2*darkkz2tilde*dsourcetilde))/darkdenominator
		darkC2=(2*exp(-darkkz1tilde*dsourcetilde)*(chi+1))/darkdenominator
		!darkA=darkC2 + darkD*exp(-2*darkkz2tilde*dsourcetilde) - exp(-darkkz1tilde*dsourcetilde)
		darkA=0 !should be zero but floating point is causing problems

		!print*, "darkD=", darkD, " darkT=",darkT," darkC2=",darkC2," darkA=", darkA
		RETURN
	END SUBROUTINE

END PROGRAM

!