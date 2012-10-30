program gaussian
	double precision kxt,kzt, thetaI, sigmat, xt, zt, pi
	double precision kxtp, dkxtp, etalim1p,etalim2p, kztp
	double precision kzt2
	double precision eta, dkxt, temp, etalim1,etalim2, eps
	double precision etatest
complex*16 field, imag, field2
complex*16 trans, ref, chi
	pi=dacos(-1.0d0)
	imag=dcmplx(0.0d0, 1.0d0)
	eta=2.0*pi*0.5
	sigmat=4.0
	thetaI=60.0d0*pi/180.0d0
	etalim1p=-eta 
	etalim2p=eta
	dkxtp=etalim2p/1000.0d0
	eps=10.0d0 

!sigmat determines beam waist

! the number 0.5=d/lamda where d is the unit of length chosen

open(20,file='absfinalhacc60cor.dat', blank='null', status='unknown')
open(21,file='absfinalhaccp60cor.dat', blank='null', status='unknown')

do xt=-20.0,20.0,0.05
	do zt=-20.0,20.0,0.05
		if(zt.lt.1.0) then
! area before interface   
field=0.0d0

do kxtp=etalim1p, etalim2p, dkxtp

	dkxt=cos(thetaI)*dkxtp
	kztp=dsqrt(eta**2-kxtp**2)
	kxt=kxtp*cos(thetaI)+kztp*sin(thetaI)
	etatest=eta**2-kxt**2
	if(etatest.ge.0.0d0) then
		kzt=dsqrt(eta**2-kxt**2)
		kzt2=dsqrt(eps*eta**2-kxt**2)

		temp=kxt*cos(thetaI)-kzt*sin(thetaI) 

		chi=dcmplx(kzt2/kzt,0.0d0)


		ref=((1.0-chi)/(1.0+chi))*cdexp(2.0*imag*kzt)

		field=field+ &
		ref*cdexp(imag*(kxt*xt-kzt*zt))*exp(-sigmat*temp**2/2.0d0)*dkxt/cos(thetaI)+ &
       cdexp(imag*(kxt*xt+kzt*zt))*exp(-sigmat*temp**2/2.0d0)*dkxt/cos(thetaI)

		endif
	end do

	else
! area after interface 
field2=0.0d0

do kxtp=etalim1p, etalim2p, dkxtp
	dkxt=cos(thetaI)*dkxtp
	kztp=dsqrt(eta**2-kxtp**2)
	kxt=kxtp*cos(thetaI)+kztp*sin(thetaI)
	etatest=eta**2-kxt**2
	if(etatest.ge.0.0d0) then
		kzt=dsqrt(eta**2-kxt**2)
		kzt2=dsqrt(eps*eta**2-kxt**2) 

		temp=kxt*cos(thetaI)-kzt*sin(thetaI) 

		chi=dcmplx(kzt2/kzt,0.0d0)


		trans=2.0d0/(1.0+chi)*cdexp(imag*(kzt-kzt2))

		field2=field2+ &
      cdexp(imag*(kxt*xt+kzt2*zt))*exp(-sigmat*temp**2/2.0d0)*dkxt/cos(thetaI)
		endif
	end do
end if 
if(zt.lt.1.0) then
	write(20,*) xt, zt, cdabs(field)
	else
		write(20,*) xt, zt, cdabs(field2)
		write(21,*) xt, zt, cdabs(field2)
		endif

	end do
end do
end program gaussian