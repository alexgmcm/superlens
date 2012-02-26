data = load('0.00002micmetal3fieldmap.dat');
i=1;
j=1;
jcount=0;


d=2E-6;
zi=-5*d;
xi=-5*d;
zf=5*d;
xf=5*d;
zstepfrac=0.1;
xstepfrac=0.1;
size=int8(((zf-zi)/(zstepfrac*d)+1));

eyarray=zeros(size,size);
eyrparray=zeros(size,size);
xarray=[xi:(xstepfrac*d):xf];
zarray=[zi:(zstepfrac*d):zf];


while (i <= size)
	while (j<=size)
		eyarray(i,j)=data((jcount* double(size) + j),3);
        eyrparray(i,j)=data((jcount* double(size) + j),4);
		j=j+1;
	end
	 jcount=jcount+1;
	 i=i+1;
	 j=1;
end


% h=surf(xarray,zarray,eyarray);
% xlabel('x');
% ylabel('z');
% zlabel('Ey');
% title('thetai=0, d=1, mu1=1, mu2=1, eps1=1, eps2=10, eta=1');
% colorbar;
% 
% %set(gca,'DataAspectRatio',[max(xarray) max(zarray) max(max(eyarray))]);
% %set(gca,'PlotBoxAspectRatio',[max(xarray) max(zarray) max(max(eyarray))]);
% 
% 
% print('-dpng', '0inc_eta1_eps11_eps210.png');

k=imagesc(xarray,zarray,eyarray);
xlabel('x');
ylabel('z');
%zlabel('Ey');
title('thetai=0pi/12, d=2microns, mu1=1, mu2=1, eps1=1, eps2=-4+2.4i, eta=1');
colorbar;
print('-dpng', 'im_0pi12inc_eta1_eps11_epsmet3c.png');

k=imagesc(xarray,zarray,eyrparray);
xlabel('x');
ylabel('z');
%zlabel('Ey');
title('real part: thetai=0pi/12, d=2microns, mu1=1, mu2=1, eps1=1, eps2=-4+2.4i, eta=1');
colorbar;
print('-dpng', 'real_im_0pi12inc_eta1_eps11_epsmet3c.png');
