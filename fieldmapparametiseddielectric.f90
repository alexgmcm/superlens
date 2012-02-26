! fieldmap.f90
! Written: 22/01/2012
! By: Alex mcMurray and Alun Daley
! Takes equations derived from EM boundary conditions
! for double interface and uses them to plot field map of Ey
! over x and z.
program fieldmap
implicit none !doesn't assume types from names, must be declared explicitly

! Declare stuff here
double precision :: x = 0.0, z = 0.0, thetai, eta, xstepfrac, zstepfrac
double precision :: xi, zi, c=300000000, xf, zf, Eyout, PI, Eyout2, eps1, eps2, mu1, mu2
complex*16 :: Ey, Ie, Re, Ce, De, Te, kx, kz1, kz2, A, B, C2, D2, F, G, H, I2, J, K, i, n1, n2, test
integer :: m, n, p !, SIZE
character :: filename*80
PI=4.D0*DATAN(1.D0) ! ensures maximum precision on any architechture apparently
i = (0.0,1.0)
!CSIN is complex?? etc.
!Can just use normal functions as modern fortran can determine the type required, 
!specialist csqrt etc. are obsolete

eta = 1 !go from 0.1 to 5, dimensionless parameter equal to d/lambda
! where d is the thickness of the slab and lambda is the free space wavelength of the incident light
! w and d and lambda are all replaced by eta

eps1=1
mu1=1
eps2=10
mu2=1



n1=SQRT(eps1*mu1)
n2=SQRT(eps2*mu2)

if (RealPart(n2) < 0) then
	n2 = -1 * n2
end if

print *, "n2real =", RealPart(n2)
print *, "n2 =", n2




!Ey=(0,0); I=(0,0); R=(0,0); C=(0,0); D=(0,0); T=(0,0); kx=(0,0); kz=(0,0)


!open(unit=2,file='fieldmap.dat') !makes output file, in unit 1

!m and n are dimensionless integers defined as integer steps of a fraction of d

!might want to use arrays as that will make it easier to keep track of values
!real, dimension(0:9) :: examplearray to make array from 0 index not 1.

!initial values of x and z
!x and z are now parametised forms equivalent to normal x and z divided by lambda
zi=-5
xi=-5

zf=5
xf=5

zstepfrac=0.1
xstepfrac=0.1
!SIZE= (((zf-zi)/zstepfrac)*((xf-xi)/xstepfrac))
!double precision, dimension(0:SIZE) :: xarray, yarray, Eyarray

do p=0, 6
	thetai = p*(PI/12.0)
	kx=n1*(2*PI*eta)*SIN(thetai)

	!kx is now a parametised kx where w/c is replaced by 2*PI*eta, 
	!therefore is multiplied by factor of d over normal kx

	kz2 = SQRT((n2*(2*PI*eta))**2 - kx**2)
	kz1 = SQRT((n1*(2*PI*eta))**2 - kx**2)

	print *, "kz1 = ", kz1


	! Remove d term from exponential arguments because this is factored into new parametised k
	! replace c/w with 1/(2*PI*eta)
	A = EXP(i*kz2); B=EXP(-i*kz2); C2=EXP(i*kz1); D2=(-1*kz1)/(2*PI*eta*mu1); F=(1*kz1)/(2*PI*eta*mu1)
	G = (-1*kz2)/(2*PI*eta*mu2); H = (1*kz2)/(2*PI*eta*mu2); I2=(-1*kz2*EXP(i*kz2))/(2*PI*eta*mu2)
	J = (1*kz2*EXP(-i*kz2))/(2*PI*eta*mu2)
	K=(-1*kz1*EXP(i*kz1))/(2*PI*eta*mu1)
	print *, "D2 = ", D2
	print *, "F = ", F
	test = (2*PI*eta*mu1)
	print *, "test = ", test 
	Ie = (1.0,0.0)
	De = ((D2-F)*(C2*I2 - A*K))/((H-F)*(A*K - C2*I2) + (G-F)*(K*B - C2*J) )
	Te = ((D2-F)*(I2*B - J*A))/((H-F)*(A*K - C2*I2) + (G-F)*(K*B - C2*J) )
	Ce = ((D2-F)*(K*B - J*C2))/((H-F)*(A*K - C2*I2) + (G-F)*(K*B - C2*J) )
	Re = Ce + De - 1.0

	print *, "Te =", Te


	write(filename,20) eta,'eta', thetai, 'thetaipdielecfieldmap.dat'
	print *, "thetai = ", thetai
	print *, "kx = ", kx 
	open(unit=2,file= filename)
	

	do m=0, int((zf-zi)/(zstepfrac))+1
		z= zi + (m*zstepfrac) !define z in terms of m and d
		!primes are dimensionless parameters
		do n=0, int(((xf-xi)/(xstepfrac)))+1
			!define x and z in their loops so they update each time
			!define x in terms of m and d
			!divide exponential terms by eta in order to compensate for parametisation
			x = xi + (n*xstepfrac)
			if (z<=0) then 
				Ey = Ie*EXP((i*kx*x)/eta)*EXP((i*kz1*z)/eta) + Re*EXP((i*kx*x)/eta)*EXP((-i*kz1*z)/eta)
				!print *, "EY: ", Ey
			else if((z>float(0)) .and. (z<eta)) then
				Ey = Ce*EXP((i*kx*x)/eta)*EXP((i*kz2*z)/eta) + De*EXP((i*kx*x)/eta)*EXP((-i*kz2*z)/eta)
			else
				Ey = Te*EXP((i*kx*x)/eta)*EXP((i*kz1*z)/eta)
			end if
			!define x and z and Ey arrays here, and fill them.
			!xarray(n+ (m*((5-xi)/xstepfrac))) = x
			!zarray(n+ (m*((5-xi)/xstepfrac))) = z
			!Eyarray(n+ (m*((5-xi)/xstepfrac))) = Ey
			Eyout = abs(Ey)
			Eyout2=RealPart(Ey)
			write(2,10) x, z, Eyout, Eyout2


		end do
			
	end do

end do
10	format(4e15.5)
20 	format(f3.1,A,f3.1,A)
print *, "kx:", kx
print *, test
end program fieldmap