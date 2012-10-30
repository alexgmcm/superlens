program newgaussfield
implicit none

double precision :: x, z, PI, sigma, kxstepfrac, n1,n2
double precision :: xi, zi, xf, zf, eta, xstepfrac, zstepfrac, phi, Ls !Ls is distance between source and interface
complex*16 :: i, kx
integer*4 :: m,n,p, xsize, zsize

double precision :: eps1, mu1, eps2, mu2
complex*16 :: Er, Et, kznorm1, kznorm2, chi

double precision, dimension(:), allocatable :: xarray
double precision, dimension(:), allocatable :: zarray

complex*16, dimension(:,:), allocatable :: fieldtransformed

complex*16, dimension(0:500) :: kxarray

xi=-20
zi=-20
xf=20
zf=20
xstepfrac = 0.1
zstepfrac = 0.1

!Ls is distance between source and interface and is what we use to parametise, setting it equal to 1
Ls=1.0d0

xsize = anint((dble(xf-xi)/xstepfrac))
zsize = anint((dble(zf-zi)/zstepfrac))

allocate(xarray(0:xsize))
allocate(zarray(0:zsize))
allocate(fieldtransformed(0:xsize,0:zsize))

PI=4.D0*DATAN(1.D0)
i=dcmplx(0.0d0,1.0d0)

phi=PI/4.0d0
eta=PI
!eta=w/c*Ls
sigma=4.0d0

!kxstepfrac=eta/100.0d0 !this should be with cos(phi)???
kxstepfrac=eta*cos(phi)/100.0d0

eps1 = 1.0d0
eps2 = 10.0d0
mu1 = 1.0d0
mu2= 1.0d0 !get NaN for normal incidence when n1=n2, look into this after

n1 = SQRT(eps1*mu1)
n2 = SQRT(eps2*mu2)

!if((real(n2)<0.and.(eps2)>0).or.&
! (eps2<0.and.(mu2<0).and.real(n2)>0))then    PUT THESE BACK IN FOR NEGATIVE REFRACTION
!  n2=-1*n2
!end if
print *, "n2=",n2
print *, "n1=",n1
 
do p=0,zsize
  zarray(p) = zi + dble(p)*zstepfrac
end do

do p=0,xsize
  xarray(p) = xi + dble(p)*xstepfrac
end do

do n=0,xsize
  do m=0, zsize
    fieldtransformed(n,m)=0
  end do
end do

do p=0,200

  kxarray(p)=-((eta)*cos(phi)) + ((dble(p)*kxstepfrac))
  kx = kxarray(p)
  kznorm1=SQRT(((n1*eta)**2.0d0) - (kx**2.0d0))
  kznorm2=SQRT(((n2*eta)**2.0d0) - (kx**2.0d0))
  
  !print *, "kznorm1=",kznorm1
  !print *, "kznorm2=",kznorm2
  
 ! if ((real(kznorm2)>0).and.(real(n2)<0)) then !this is right   PUT THIS BACK IN AFTER
 !  kznorm2 = -kznorm2
 ! end if
	
 ! if (real(kznorm1)<0) then   PUT THIS BACK IN
 !  kznorm1 = -kznorm1
 ! end if

   chi=(kznorm2*mu1)/(kznorm1*mu2)

   Er = (cdEXP(2*i*kznorm1*Ls))*(1-chi)/(1+chi)
   Et = (cdEXP(i*(kznorm1-kznorm2)*Ls))*2.0d0/(1+chi)

  open(unit=2,file='gfield020.dat')

  do m=0, zsize
    z=zarray(m)
  
    do n=0, xsize
      x=xarray(n)
      kx=kxarray(p)

  !zrofra = (zarray(m)*DCOS(phi))+(xarray(n)*DSIN(phi))
  !xrofra = (xarray(n)*DCOS(phi))-(zarray(m)*DSIN(phi))
  
  if(z<=Ls) then !this is where we say the source is at z=0, defined negatively in else
    fieldtransformed(n,m)=fieldtransformed(n,m)+ &
     ((dSQRT(sigma/(2.0d0*PI))*cdEXP(i*kx*x)*cdEXP(i*kznorm1*z)*cdEXP(-(sigma/2.0d0)* &
      (((kx*cos(phi))-(kznorm1*sin(phi)))**2.0d0))*kxstepfrac/cos(phi))+&
       (Er*dSQRT(sigma/(2.0d0*PI))*cdEXP(i*kx*x)*cdEXP(-i*kznorm1*z)*cdEXP(-(sigma/2.0d0)* &
      (((kx*cos(phi))-(kznorm1*sin(phi)))**2.0d0))*kxstepfrac/cos(phi)))
  else
    fieldtransformed(n,m)=fieldtransformed(n,m)+ &
     (Et*dSQRT(sigma/(2.0d0*PI))*cdEXP(i*kx*x)*cdEXP(i*kznorm2*z)*cdEXP(-(sigma/2.0d0)* &
      (((kx*cos(phi))-(kznorm1*sin(phi)))**2.0d0))*kxstepfrac/cos(phi))!change back to 2
  end if
!could be a problem with the numerical integration method, which is why reflection isn't shown. Put it back to normal expressions(factorise)

    end do
  end do
end do

do m=0,zsize
  do n=0,xsize
  
  write(2,*) xarray(n),zarray(m), abs(fieldtransformed(n,m)), real(fieldtransformed(n,m))
 
  end do
end do

end program newgaussfield