
zi=-20;
xi=-20;
zf=20;
xf=20;
zstepfrac=0.1;
xstepfrac=0.1;
size=ceil(((zf-zi)/(zstepfrac)))+1;
sizex=ceil(((xf-xi)/(xstepfrac)))+1;
eta='pi';
eps2=-1.0;
mu2=-1.0;
thetai='60';
dsource=1;
 
g=[1.0];
for x=1:1
g(x) = round(g(x)*10^1)/(10^1);
end
sizesquare=ceil(double(size)*double(sizex));
data=zeros(sizesquare, 4, 7);

for x=1:1
num2str(g(x), '%3.1f');
data(:,:,x) = load(strcat('data/singint',thetai,'degs','3.1eta4.0sigmatilde.dat'));
end

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
xlabel('x/lambda');
ylabel('z/lambda');
title(strcat('g=',num2str(g(1)),',thetai=',num2str(thetai),',mu1=1, mu2=',num2str(mu2),', eps1=1, eps2=',num2str(eps2),', eta=',num2str(eta) ));
colorbar;
line([xi xf],[dsource dsource],'linewidth',4,'Color', 'k');
%line([xi xf],[0 0],'linewidth',4,'Color', 'k');
print('-dpng',strcat('plots/im_singint',thetai,'degspieta4.0sigmatilde',num2str(eps2),'eps2.png'));


k=imagesc(xarray,zarray,eyrparray);
xlabel('x/lambda');
ylabel('z/lambda');
title(strcat('real part: g=',num2str(g(1)),',thetai=',num2str(thetai),',mu1=1, mu2=',num2str(mu2),', eps1=1, eps2=',num2str(eps2),', eta=',num2str(eta) ));
colorbar;
line([xi xf],[dsource dsource],'linewidth',4,'Color', 'k');
%line([xi xf],[0 0],'linewidth',4,'Color', 'k');
print('-dpng',strcat('plots/real_singint',thetai,'degspieta4.0sigmatilde',num2str(eps2),'eps2.png'));
end
