
zi=-5;
xi=-5;
zf=5;
xf=5;
zstepfrac=0.1;
xstepfrac=0.1;
size=int8(((zf-zi)/(zstepfrac)+1));
eta = 1;
eps2=10;

 
angle=[0.0:(pi/12.0):pi/2.0]; 
for x=1:7
angle(x) = round(angle(x)*10^1)/(10^1);
end
sizesquare=int32(double(size)*double(size));
data=zeros(sizesquare, 4, 7);

for x=1:7
num2str(angle(x), '%3.1f');
data(:,:,x) = load(strcat(num2str(eta, '%3.1f'),'eta',num2str(angle(x), '%3.1f'),'thetaipdielecfieldmap.dat'));
end

for x=1:7

i=1;
j=1;
jcount=0;



eyarray=zeros(size,size);
eyrparray=zeros(size,size);
xarray=[xi:(xstepfrac):xf];
zarray=[zi:(zstepfrac):zf];


while (i <= size)
	while (j<=size)
		eyarray(i,j)=data((jcount* double(size) + j),3,x);
        eyrparray(i,j)=data((jcount* double(size) + j),4,x);
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
xlabel('x/lambda');
ylabel('z/lambda');
%zlabel('Ey');
title(strcat('thetai=',num2str((x-1)),'pi/12, mu1=1, mu2=1, eps1=1, eps2=',num2str(eps2),', eta=',num2str(eta) ));
colorbar;
line([xi xf],[eta eta],'linewidth',4,'Color', 'k');
line([xi xf],[0 0],'linewidth',4,'Color', 'k');
print('-dpng', strcat('im_',num2str((x-1)),'pi12inc_eta',num2str(eta),'_eps11_eps2',num2str(eps2),'.png'));


k=imagesc(xarray,zarray,eyrparray);
xlabel('x/lambda');
ylabel('z/lambda');
%zlabel('Ey');
title(strcat('real part: thetai=',num2str((x-1)),'pi/12, mu1=1, mu2=1, eps1=1, eps2=',num2str(eps2),', eta=',num2str(eta) ));
colorbar;
line([xi xf],[eta eta],'linewidth',4,'Color', 'k');
line([xi xf],[0 0],'linewidth',4,'Color', 'k');
print('-dpng', strcat('real_im_',num2str((x-1)),'pi12inc_eta',num2str(eta),'_eps11_eps2',num2str(eps2),'.png'));
end