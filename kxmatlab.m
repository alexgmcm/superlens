kxdata(:,:)=load('data/kxcombinedimageline90.0degs3.1eta0.010sigmatilde3.0secint30etalimit.dat');
kx=kxdata(:,1);
r=kxdata(:,2);
c=kxdata(:,3);
d=kxdata(:,4);
t=kxdata(:,5);

plot(kx,r,'xr');
hold on;
plot(kx,c,'xb');
plot(kx,d,'.g');
plot(kx,t,'xm');
xlabel('kx');
ylabel('coefficients');
legend('r','c','d','t');
print('-dpng','plots/kxcoeffsigma0.010.png');
