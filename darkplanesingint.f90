program darkplanesingint
implicit none !don't assume types
double precision :: x = 0.0, z = 0.0, thetai, eta, xstepfrac, zstepfrac
double precision :: xi, zi, c=300000000, xf, zf, Eyout, PI, Eyout2, eps1, eps2, mu1, mu2, kx
complex*16 :: Ey, Re, Te, kz1, kz2, i, n1, n2, test, chi
integer :: m, n, p !, SIZE
character :: filename*80

PI=4.D0*DATAN(1.D0) ! ensures maximum precision on any architechture apparently
i = (0.0,1.0)
eta = PI

eps1=1
mu1=1
eps2=10
mu2=1



n1=SQRT(eps1*mu1)
n2=SQRT(eps2*mu2)

if (RealPart(n2) < 0) then
	n2 = -1 * n2
end if


zi=-5
xi=-5

zf=5
xf=5

zstepfrac=0.1
xstepfrac=0.1 

kx=1.5

kz2 = SQRT((n2*(eta))**2 - kx**2)
kz1 = SQRT((n1*(eta))**2 - kx**2)

chi=(kz2*mu1)/(kz1*mu2)


Re = exp(-kz1)*(2/(chi+1))
Te = exp(-kz1)*((2/(chi+1)) - 1)

write(filename,20) 'data/darksingint', kx,'kx.dat'
open(unit=2,file= filename)


	do m=0, nint(((zf-zi)/zstepfrac))
		z= zi + (m*zstepfrac) !define z in terms of m and d
		!primes are dimensionless parameters
		do n=0, nint(((xf-xi)/xstepfrac))
			!define x and z in their loops so they update each time
			!define x in terms of m and d
			!divide exponential terms by eta in order to compensate for parametisation
			x = xi + (n*xstepfrac)
			if (z<=1) then 
				Ey = EXP((i*kx*x))*EXP((-kz1*z)) + Re*EXP((i*kx*x))*EXP((kz1*(z - 1)))
				!print *, "EY: ", Ey
			else
				Ey = Te*EXP((i*kx*x))*EXP((-kz2*(z - 1)))
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


10	format(4e15.5,4e15.5,4e15.5,4e15.5)
20 	format(A,f3.1,A)





end program darkplanesingint