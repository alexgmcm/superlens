function resolutionsuperlenscliscript(secondinterface, thetamax)
%eyabs is now intensity
zi=-20;
xi=-20;
zf=20;
xf=20;
zstepfrac=0.1;
xstepfrac=0.1;
size=ceil(((zf-zi)/(zstepfrac)))+1;
sizex=ceil(((xf-xi)/(xstepfrac)))+1;
eta='pi';
eps2=-1;
mu2=-1;
thetai='0';
dsource=1;
%secondinterface=3*dsource;
 
g=0.01
%for x=1:1
%g(x) = round(g(x)*10^1)/(10^1);
%end
sizesquare=ceil(double(size)*double(sizex));
data=zeros(sizesquare, 4, 7);


gs=num2str(g, '%3.1f');
data(:,:,1) = load(strcat('data/res',num2str(thetamax),'degs','3.1eta',gs,'sigmatilde',num2str(secondinterface),'secint.dat'
open(unit=2,file= filename).dat'));


for x=1:1

i=1;
j=1;
jcount=0;



eyarray=zeros(sizex,size);%this is the tranformed field modulus
eyrparray=zeros(sizex,size);%this is the transformed field real part
xarray=[xi:(xstepfrac):xf];
xzeroindex=find(xarray==0);
zarray=[zi:(zstepfrac):zf];
zzeroindex=find(zarray==0);
zsecondinterface=find(zarray==secondinterface)


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
title(strcat('g=',gs,',thetai=',num2str(thetai),'degs,mu1=1, mu2=',num2str(mu2),', eps1=1, eps2=',num2str(eps2),', eta=',eta ));
colorbar;
line([xi xf],[dsource dsource],'linewidth',4,'Color', 'k');
line([xi xf],[secondinterface secondinterface],'linewidth',4,'Color', 'k');
%line([xi xf],[0 0],'linewidth',4,'Color', 'k');
print('-dpng',strcat('plots/im_res',num2str(thetamax),'degspieta',gs,'sigmatilde',num2str(eps2),'eps2',num2str(secondinterface),'secint.png'));


k=imagesc(xarray,zarray,eyrparray);
xlabel('x/dsource');
ylabel('z/dsource');
title(strcat('real part: g=',gs,',thetai=',num2str(thetai),'degs,mu1=1, mu2=',num2str(mu2),', eps1=1, eps2=',num2str(eps2),', eta=',eta ));
colorbar;
line([xi xf],[dsource dsource],'linewidth',4,'Color', 'k');
line([xi xf],[secondinterface secondinterface],'linewidth',4,'Color', 'k');
%line([xi xf],[0 0],'linewidth',4,'Color', 'k');
print('-dpng',strcat('plots/real_res',num2str(thetamax),'degspieta',gs,'sigmatilde',num2str(eps2),'eps2',num2str(secondinterface),'secint.png'));


imageaxis=eyarray(xzeroindex,:);
imageaxis=imageaxis(zsecondinterface:length(imageaxis));
[maxval,imageindex]=max(imageaxis);
zimagepos=zarray(imageindex);

slabthickness=secondinterface-dsource;


ximageaxis=eyarray(:,zimagepos)
[fhwmindices]=find((ximageaxis==(max(imageaxis)/2)),2);
fwhm=xarray(fwhmindices(2))-xarray(fwhmindices(1))









end
