! fieldmap.f90
! Written: 22/01/2012
! By: Alex mcMurray and Alun Daley
! Takes equations derived from EM boundary conditions
! for double interface and uses them to plot field map of Ey
! over x and z.




!complex function getefield(x,z) result(Ey)
!
	
!end function getefield


program fieldmap
implicit none !doesn't assume types from names, must be declared explicitly

! Declare stuff here
double precision :: x = 0.0, z = 0.0, d = 0.0, thetai, eta, xstepfrac, zstepfrac
double precision :: xi, zi, c=300000000, xf, zf, w, Eyout, PI, Eyout2, test, eps1, eps2, mu1, mu2
complex*16 :: Ey, Ie, Re, Ce, De, Te, kx, kz1, kz2, A, B, C2, D2, F, G, H, I2, J, K, i, n1, n2
integer :: m, n, p !, SIZE
character :: filename*80
PI=4.D0*DATAN(1.D0) ! ensures maximum precision on any architechture apparently
i = (0.0,1.0)
!CSIN is complex?? etc.
!Can just use normal functions as modern fortran can determine the type required, 
!specialist csqrt etc. are obsolete
d=2E-6 !choose value
eta = 1 !go from 0.1 to 5
w = c*eta/d !eta is dimensionless parameter 
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


kx=n1*(w/c)*SIN(thetai)
kz2 = SQRT((n2*w/c)**2 - kx**2)
kz1 = SQRT((n1*w/c)**2 - kx**2)


! ELIMINATE w for DIMENSIONLESS VARIABLE
A = EXP(i*kz2*d); B=EXP(-i*kz2*d); C2=EXP(i*kz1*d); D2=(-c*kz1)/(w*mu1); F=(c*kz1)/(w*mu1)
G = (-c*kz2)/(w*mu2); H = (c*kz2)/(w*mu2); I2=(-c*kz2*EXP(i*kz2*d))/(w*mu2)
J = (c*kz2*EXP(-i*kz2*d))/(w*mu2)
K=(-c*kz1*EXP(i*kz1*d))/(w*mu1)


Ie = (1.0,0.0)
De = ((D2-F)*(C2*I2 - A*K))/((H-F)*(A*K - C2*I2) + (G-F)*(K*B - C2*J) )
Te = ((D2-F)*(I2*B - J*A))/((H-F)*(A*K - C2*I2) + (G-F)*(K*B - C2*J) )
Ce = ((D2-F)*(K*B - J*C2))/((H-F)*(A*K - C2*I2) + (G-F)*(K*B - C2*J) )
Re = Ce + De - 1.0

print *, "Te =", Te



!Ey=(0,0); I=(0,0); R=(0,0); C=(0,0); D=(0,0); T=(0,0); kx=(0,0); kz=(0,0)


!open(unit=2,file='fieldmap.dat') !makes output file, in unit 1

!m and n are dimensionless integers defined as integer steps of a fraction of d

!might want to use arrays as that will make it easier to keep track of values
!real, dimension(0:9) :: examplearray to make array from 0 index not 1.

!initial values of x and z
zi=-5*d
xi=-5*d

zf=5*d
xf=5*d

zstepfrac=0.1
xstepfrac=0.1
!SIZE= (((zf-zi)/zstepfrac)*((xf-xi)/xstepfrac))
!double precision, dimension(0:SIZE) :: xarray, yarray, Eyarray

do p=0, 6
	thetai = p*PI/12.0
	write(filename,20) thetai, 'dielecfieldmap.dat'
	print *, thetai
	open(unit=2,file= filename)
	

	do m=0, int((zf-zi)/(zstepfrac*d))+1
		z= zi + (m*zstepfrac*d) !define z in terms of m and d
		!primes are dimensionless parameters
		do n=0, int(((xf-xi)/(xstepfrac*d)))+1
			!define x and z in their loops so they update each time
			!define x in terms of m and d
			x = xi + (n*xstepfrac * d)
			if (z<=0) then 
				Ey = Ie*EXP((i*kx*x))*EXP(i*kz1*z) + Re*EXP(i*kx*x)*EXP(-i*kz1*z)
				!print *, "EY: ", Ey
			else if((z>float(0)) .and. (z<d)) then
				Ey = Ce*EXP(i*kx*x)*EXP(i*kz2*z) + De*EXP(i*kx*x)*EXP(-i*kz2*z)
			else
				Ey = Te*EXP(i*kx*x)*EXP(i*kz1*z)
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
20 	format(f6.4,A)

print *, "kx:", kx
print *, test
end program fieldmap


!define ey data as real, modulus of ey then take into array. check mathematics between types to generate ey
! fix upper limits on loops











