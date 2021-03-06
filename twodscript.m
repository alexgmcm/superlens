function twodscript(secondinterface, thetamax, g, string)
%eyabs is now intensity
zi=0;
xi=-10;
zf=10;
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
kxcutoff='30';
%secondinterface=3*dsource;
firstinterface=dsource;
imageplane=4*dsource;

%g=0.001
%for x=1:1
%g(x) = round(g(x)*10^1)/(10^1);
%end
sizesquare=ceil(double(size)*double(sizex));
data=zeros(sizesquare, 4, 7);


gs=num2str(g, '%5.3f');
data(:,:,1) = load(strcat('data/',string,'2D',num2str(thetamax, '%3.1f'),'degs','3.1eta',gs,'sigmatilde',num2str(secondinterface,'%2.1f'),'secint30etalimit.dat'));



for x=1:1

i=1;
j=1;
jcount=0;



intensarray=zeros(size,sizex);%this is the tranformed field modulus
eyarray=zeros(size,sizex);%this is the transformed field real part
xarray=[xi:(xstepfrac):xf];
xzeroindex=find(xarray==0);
zarray=[zi:(zstepfrac):zf];
zzeroindex=find(zarray==0);
zsecondinterface=find(zarray==secondinterface);
zfirstinterface=find(zarray==firstinterface);


while (i <= size)
	while (j<=sizex)
		intensarray(i,j)=data((jcount* double(sizex) + j),3,x);
        eyarray(i,j)=data((jcount* double(sizex) + j),4,x);
		j=j+1;
	end
	 jcount=jcount+1;
	 i=i+1;
	 j=1;
end

k=imagesc(xarray,zarray,log(intensarray)); %log or cut off maxcolor???
xlabel('x/dsource','FontSize',16);
ylabel('z/dsource','FontSize',16);
%title(strcat('g=',gs,',thetai=',num2str(thetai),'degs,mu1=1, mu2=',num2str(mu2),', eps1=1, eps2=',num2str(eps2),', eta=',eta,' kxcutoff=',kxcutoff ));
title(strcat('logarithmic plot of intensity, sigma=',gs,',  ',string ),'FontSize',16);
colorbar;
caxis([-10 20]);
line([xi xf],[dsource dsource],'linewidth',4,'Color', 'k');
line([xi xf],[secondinterface secondinterface],'linewidth',4,'Color', 'k');
%line([xi xf],[0 0],'linewidth',4,'Color', 'k');
print('-depsc2',strcat('plots/',string,'2D',num2str(thetamax),'degspieta',gs,'sigmatilde',num2str(eps2),'eps2','kxcutoff',kxcutoff,num2str(secondinterface),'secint.eps'));


zimageaxis=eyarray(:,xzeroindex);
zimagepos=find(zarray==imageplane);


% plot(zarray,zimageaxis,'-r');
% title(strcat('Intensity vs z, sigma=',gs,' ',string));
% xlabel('z ,d_s');
% ylabel('Intensity, EE*');
% line([dsource dsource], [0 ceil(max(zimageaxis))],'linewidth',1,'Color', 'k');
% line([secondinterface secondinterface],[0 ceil(max(zimageaxis))],'linewidth',1,'Color', 'k');
% line([zimagepos zimagepos], [0 ceil(max(zimageaxis))],'linewidth',1,'Color', 'k');
% line([0 max(zarray)], [zimageaxis(1) zimageaxis(1)],'linewidth',1,'Color', 'k','LineStyle','--');
% ylim([0 10]);
% xlim([0 10]);
% print('-dpng',strcat('plots/',string,'eyprofilegs',gs,'kxcutoff',kxcutoff,'.png'));


end