function resolutionsuperlenscliscript(secondinterface, thetamax)
%eyabs is now intensity
zi=-5;
xi=-10;
zf=15;
xf=10;
zstepfrac=0.1;
xstepfrac=0.01;
size=ceil(((zf-zi)/(zstepfrac)))+1;
sizex=ceil(((xf-xi)/(xstepfrac)))+1;
eta='pi';
eps2=-1;
mu2=-1;
thetai='0';
dsource=1;
%secondinterface=3*dsource;
 
g=0.001
%for x=1:1
%g(x) = round(g(x)*10^1)/(10^1);
%end
sizesquare=ceil(double(size)*double(sizex));
data=zeros(sizesquare, 4, 7);


gs=num2str(g, '%5.3f');
data(:,:,1) = load(strcat('data/res',num2str(thetamax, '%3.1f'),'degs','3.1eta',gs,'sigmatilde',num2str(secondinterface,'%2.1f'),'secint.dat'));



for x=1:1

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


imageaxis=eyarray(:,xzeroindex);
imageaxis=imageaxis(zsecondinterface:length(imageaxis));
[maxval,imageindex]=max(imageaxis);
zimagepos=zarray(zsecondinterface+imageindex-1);

slabthickness=secondinterface-dsource;

ximageaxis=eyarray(zsecondinterface+imageindex,:);
xsourceaxis=eyarray(zzeroindex,:);

plot(xarray,ximageaxis,'-r');
title(strcat('Intensity profile (EE*) at image, thetamax: ', num2str(thetamax),' degs'));
xlabel('x ,d_s');
ylabel('Intensity (EE*)');
%ylim([0 7e-03])
print('-dpng',strcat('plots/imageprofile',num2str(thetamax),'degsthetamax.png'));


plot(xarray,xsourceaxis,'-r');
title(strcat('Intensity profile (EE*) at source, thetamax: ', num2str(thetamax),' degs'));
xlabel('x ,d_s');
ylabel('Intensity (EE*)');
print('-dpng',strcat('plots/sourceprofile',num2str(thetamax),'degsthetamax.png'));

[xmaxval,ximageindex]=max(ximageaxis);
lhs=ximageaxis(1:ximageindex);
rhs=ximageaxis(ximageindex:length(ximageaxis));
llowerind= find(lhs<=(max(imageaxis)/2),1,'last');
lhigherind=find(lhs>=(max(imageaxis)/2),1,'first');
leftval=(xarray(lhigherind)+xarray(llowerind))/2;

rlowerind=find(rhs>=(max(imageaxis)/2),1,'last');
rhigherind= find(rhs<=(max(imageaxis)/2),1,'first');
rightval=(xarray(ximageindex+rhigherind-1)+xarray(ximageindex+rlowerind-1))/2;

fwhm=rightval-leftval;
expectedpendry=(slabthickness-dsource)+secondinterface;

fid = fopen(strcat('data/plotdata',num2str(secondinterface),'secint.txt'), 'a');
outdata = [slabthickness;zimagepos;fwhm;expectedpendry;thetamax];
fprintf(fid, '%3.1f %3.1f %6.4f %3.1f %3.1f\n', outdata);
fclose(fid);







end
