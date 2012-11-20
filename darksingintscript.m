function darksingintscript(kx)

zi=0;
xi=0;
zf=3;
xf=3;
zstepfrac=0.01;
xstepfrac=0.01;
size=ceil(((zf-zi)/(zstepfrac)))+1;
sizex=ceil(((xf-xi)/(xstepfrac)))+1;
eta='pi';
eps2=10;
mu2=1;
thetai='0';
dsource=1;
%secondinterface=3*dsource;
 
%g=0.01
%for x=1:1
%g(x) = round(g(x)*10^1)/(10^1);
%end
sizesquare=ceil(double(size)*double(sizex));
data=zeros(sizesquare, 4, 7);


%gs=num2str(g, '%3.1f');
data(:,:,1) = load(strcat('data/darksingint',num2str(kx),'kx.dat'));


for x=1:1

i=1;
j=1;
jcount=0;



eyarray=zeros(sizex,size);%this is the tranformed field modulus
eyrparray=zeros(sizex,size);%this is the transformed field real part
xarray=[xi:(xstepfrac):xf];
zarray=[zi:(zstepfrac):zf];


while (i <= sizex)
	while (j<=size)
		eyarray(i,j)=data((jcount* double(size) + j),3,x);
        eyrparray(i,j)=data((jcount* double(size) + j),4,x);
		j=j+1;
	end
	 jcount=jcount+1;
	 i=i+1;
	 j=1;
end

k=imagesc(xarray,zarray,eyarray);
xlabel('x/dsource');
ylabel('z/dsource');
title(strcat(',thetai=',num2str(thetai),'degs,mu1=1, mu2=',num2str(mu2),', eps1=1, eps2=',num2str(eps2),', eta=',eta,' kx=',num2str(kx)  ));
colorbar;
line([xi xf],[dsource dsource],'linewidth',4,'Color', 'k');
%line([xi xf],[0 0],'linewidth',4,'Color', 'k');
print('-dpng',strcat('plots/im_darksingint',thetai,'degspieta',num2str(kx),'kx',num2str(eps2),'eps2.png'));


k=imagesc(xarray,zarray,eyrparray);
xlabel('x/dsource');
ylabel('z/dsource');
title(strcat('real part:,thetai=',num2str(thetai),'degs,mu1=1, mu2=',num2str(mu2),', eps1=1, eps2=',num2str(eps2),', eta=',eta,' kx=',num2str(kx) ));
colorbar;
line([xi xf],[dsource dsource],'linewidth',4,'Color', 'k');
%line([xi xf],[0 0],'linewidth',4,'Color', 'k');
print('-dpng',strcat('plots/real_darksingint',thetai,'degspieta',num2str(kx),'kx',num2str(eps2),'eps2.png'));
end
