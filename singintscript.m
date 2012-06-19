
zi=0;
xi=-10;
zf=20;
xf=10;
zstepfrac=0.1;
xstepfrac=0.1;
size=int32(((zf-zi)/(zstepfrac)+1));
eta=1;
eps2=5.0;
mu2=1.0;
thetai='PI/4.0';
dsource=10;
 
g=[1.0];
for x=1:1
g(x) = round(g(x)*10^1)/(10^1);
end
sizesquare=int32(double(size)*double(size));
data=zeros(sizesquare, 4, 7);

for x=1:1
num2str(g(x), '%3.1f');
data(:,:,x) = load('data/singint0.16pirads1.0eta1.0sigmatilde.dat');
end

for x=1:1

i=1;
j=1;
jcount=0;



eyarray=zeros(size,size);%this is the tranformed field modulus
eyrparray=zeros(size,size);%this is the transformed field real part
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

k=imagesc(xarray,zarray,eyarray);
xlabel('x/lambda');
ylabel('z/lambda');
title(strcat('g=',num2str(g(1)),',thetai=',num2str(thetai),',mu1=1, mu2=',num2str(mu2),', eps1=1, eps2=',num2str(eps2),', eta=',num2str(eta) ));
colorbar;
line([xi xf],[dsource dsource],'linewidth',4,'Color', 'k');
%line([xi xf],[0 0],'linewidth',4,'Color', 'k');
print('-dpng','plots/im_singint0.16pirads1.0eta1.0sigmatilde.png');


k=imagesc(xarray,zarray,eyrparray);
xlabel('x/lambda');
ylabel('z/lambda');
title(strcat('real part: g=',num2str(g(1)),',thetai=',num2str(thetai),',mu1=1, mu2=',num2str(mu2),', eps1=1, eps2=',num2str(eps2),', eta=',num2str(eta) ));
colorbar;
line([xi xf],[dsource dsource],'linewidth',4,'Color', 'k');
%line([xi xf],[0 0],'linewidth',4,'Color', 'k');
print('-dpng','plots/real_singint0.16pirads1.0eta1.0sigmatilde.png');
end
