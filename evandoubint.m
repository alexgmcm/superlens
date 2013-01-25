
%eyabs is now intensity
zi=0;
xi=-10;
zf=20;
xf=10;
zstepfrac=0.1;
xstepfrac=0.01;
size=ceil(((zf-zi)/(zstepfrac)))+1;
sizex=ceil(((xf-xi)/(xstepfrac)))+1;
eta='1';
kx=1.5;
eps2=-1;
mu2=-1;
thetai='0';
dsource=1;
secondinterface=3*dsource;

sizesquare=ceil(double(size)*double(sizex));
data=zeros(sizesquare, 4, 7);



data(:,:,1) = load(strcat('data/evandoubint1.5kx1.0eta3.0secint.dat'));



x=1;

i=1;
j=1;
jcount=0;



eyarray=zeros(size,sizex);%this is the tranformed field modulus
eyrparray=zeros(size,sizex);%this is the transformed field real part
xarray=[xi:(xstepfrac):xf];
xzeroindex=find(xarray==0);
zarray=[zi:(zstepfrac):zf];
zzeroindex=find(zarray==0);
zsecondinterface=find(zarray==secondinterface);


while (i <= size)
	while (j<=sizex)
		eyarray(i,j)=data((jcount* double(sizex) + j),3,x);
        eyrparray(i,j)=data((jcount* double(sizex) + j),4,x);
		j=j+1;
	end
	 jcount=jcount+1;
	 i=i+1;
	 j=1;
end

k=imagesc(xarray,zarray,eyarray);
xlabel('x/dsource');
ylabel('z/dsource');
title(strcat('kx=',num2str(kx),',mu1=1, mu2=',num2str(mu2),', eps1=1, eps2=',num2str(eps2),', eta=',eta ));
colorbar;
line([xi xf],[dsource dsource],'linewidth',4,'Color', 'k');
line([xi xf],[secondinterface secondinterface],'linewidth',4,'Color', 'k');
%line([xi xf],[0 0],'linewidth',4,'Color', 'k');
print('-dpng',strcat('plots/evandoubint',num2str(kx),'kx',eta,'eta.png'));


k=imagesc(xarray,zarray,eyrparray);
xlabel('x/dsource');
ylabel('z/dsource');
title(strcat('kx=',num2str(kx),',mu1=1, mu2=',num2str(mu2),', eps1=1, eps2=',num2str(eps2),', eta=',eta ));
colorbar;
line([xi xf],[dsource dsource],'linewidth',4,'Color', 'k');
line([xi xf],[secondinterface secondinterface],'linewidth',4,'Color', 'k');
%line([xi xf],[0 0],'linewidth',4,'Color', 'k');
print('-dpng',strcat('plots/real_evandoubint',num2str(kx),'kx',eta,'eta.png'));


imageaxis=eyarray(:,xzeroindex);


plot(zarray,imageaxis,'-r');
title(strcat('abs(Ey) vs z, kx:', num2str(kx),' eta:', eta ));
xlabel('z ,d_s');
ylabel('abs(Ey)');
%ylim([0 7e-03])
print('-dpng',strcat('plots/eyrpprofile',num2str(kx),'kxeta',eta,'.png'));

