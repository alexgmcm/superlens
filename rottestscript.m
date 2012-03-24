
zi=-10;
xi=-10;
zf=10;
xf=10;
zstepfrac=0.1;
xstepfrac=0.1;
size=int32(((zf-zi)/(zstepfrac)+1));
eta = 3;
thetai='0.25pi';

 
g=[1.0]; 
for x=1:1
g(x) = round(g(x)*10^1)/(10^1);
end
sizesquare=int32(double(size)*double(size));
data=zeros(sizesquare, 6, 7);

for x=1:1
num2str(g(x), '%3.1f');
data(:,:,x) = load(strcat('data/rottest',thetai,'rads',num2str(eta, '%3.1f'),'eta',num2str(g(x), '%3.1f'),'sigmatilde.dat'));
end

for x=1:1

i=1;
j=1;
jcount=0;



eyrealarray=zeros(size,size);%this is the tranformed field modulus
eykspacearray=zeros(size,size);%this is the transformed field real part
xarray=[xi:(xstepfrac):xf];
zarray=[zi:(zstepfrac):zf];


while (i <= size)
	while (j<=size)
		eyrealarray(i,j)=data((jcount* double(size) + j),3,x);
        eyrealarrayrp(i,j)=data((jcount* double(size) + j),5,x);
        eykspacearray(i,j)=data((jcount* double(size) + j),4,x);
        eykspacearrayrp(i,j)=data((jcount* double(size) + j),6,x);
		j=j+1;
	end
	 jcount=jcount+1;
	 i=i+1;
	 j=1;
end

k=imagesc(xarray,zarray,eyrealarray);
xlabel('x/d');
ylabel('z/d');
%zlabel('Ey');
title(strcat('real: sigmatilde=',num2str(g(1)),',thetai=',num2str(thetai),', eta=',num2str(eta) ));
colorbar;
%line([xi xf],[eta eta],'linewidth',4,'Color', 'k');
%line([xi xf],[0 0],'linewidth',4,'Color', 'k');
print('-dpng', strcat('plots/rottestreal_thetai',thetai,'_g',num2str(g(1)),'_eta',num2str(eta),'.png'));

k=imagesc(xarray,zarray,eyrealarrayrp);
xlabel('x/d');
ylabel('z/d');
%zlabel('Ey');
title(strcat('realrp: sigmatilde=',num2str(g(1)),',thetai=',num2str(thetai),', eta=',num2str(eta) ));
colorbar;
%line([xi xf],[eta eta],'linewidth',4,'Color', 'k');
%line([xi xf],[0 0],'linewidth',4,'Color', 'k');
print('-dpng', strcat('plots/rottestrealrp_thetai',thetai,'_g',num2str(g(1)),'_eta',num2str(eta),'.png'));



k=imagesc(xarray,zarray,eykspacearray);
xlabel('x/d');
ylabel('z/d');
%zlabel('Ey');
title(strcat('fourier: sigmatilde=',num2str(g(1)),',thetai=',num2str(thetai),', eta=',num2str(eta) ));
colorbar;
%line([xi xf],[eta eta],'linewidth',4,'Color', 'k');
%line([xi xf],[0 0],'linewidth',4,'Color', 'k');
print('-dpng', strcat('plots/rottestkspace_thetai',thetai,'_g',num2str(g(1)),'_eta',num2str(eta),'.png'));


k=imagesc(xarray,zarray,eykspacearrayrp);
xlabel('x/d');
ylabel('z/d');
%zlabel('Ey');
title(strcat('fourierrp: sigmatilde=',num2str(g(1)),',thetai=',num2str(thetai),', eta=',num2str(eta) ));
colorbar;
%line([xi xf],[eta eta],'linewidth',4,'Color', 'k');
%line([xi xf],[0 0],'linewidth',4,'Color', 'k');
print('-dpng', strcat('plots/rottestkspacerp_thetai',thetai,'_g',num2str(g(1)),'_eta',num2str(eta),'.png'));

end