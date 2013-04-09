function masterbetamatlab(sigma)

propdata(:,:,1) = load(strcat('data/proponlyimageline90.0degs3.1eta',num2str(sigma,'%5.3f'),'sigmatilde3.0secint30etalimit.dat'));
propx=propdata(:,1);
propintensity=propdata(:,3);
propefield=propdata(:,4);

combdata(:,:,1) = load(strcat('data/combinedimageline90.0degs3.1eta',num2str(sigma,'%5.3f'),'sigmatilde3.0secint30etalimit.dat'));
combx=combdata(:,1);
combintensity=combdata(:,3);
combefield=combdata(:,4);

darkdata(:,:,1) = load(strcat('data/darkonlyimageline90.0degs3.1eta',num2str(sigma,'%5.3f'),'sigmatilde3.0secint30etalimit.dat'));
darkx=darkdata(:,1);
darkintensity=darkdata(:,3);
darkefield=darkdata(:,4);

plot(propx,propintensity,'-r');
title(strcat('Intensity versus x (in dsource), symmetric case, sigma=',num2str(sigma,'%5.3f')));
xlabel('x / dsource');
ylabel('Intensity , EE*');
hold on;
plot(combx,combintensity,'-b');
plot(darkx,darkintensity,'-g');
legend('propagating only','combined', 'dark');
print('-dpng',strcat('plots/intensityoverlaysigma',num2str(sigma,'%5.3f'),'.png'));
hold off;

plot((propx/2),propintensity,'-r','LineWidth',4);
title(strcat('Intensity versus x (in lambda), symmetric case, sigma=',num2str(sigma,'%5.3f')),'FontSize',16);
xlabel('x / $\lambda$', 'Interpreter', 'LaTex','FontSize',16);
ylabel('Intensity , EE*','FontSize',16);
hold on;
plot((combx/2),combintensity,'-b','LineWidth',4);
plot((darkx/2),darkintensity,'-g','LineWidth',4);
leg=legend('propagating only','combined', 'dark')
set(leg,'FontSize',16)
print('-depsc2',strcat('plots/intensityoverlaylambdasigma',num2str(sigma,'%5.3f'),'.eps'));
hold off;

plot(propx,propefield,'-r');
title(strcat('E_y versus x (in dsource), symmetric case, sigma=',num2str(sigma,'%5.3f')));
xlabel('x / dsource');
ylabel('E_y');
hold on;
plot(combx,combefield,'-b');
plot(darkx,darkefield,'-g');
legend('propagating only','combined', 'dark');
print('-dpng',strcat('plots/efieldoverlaysigma',num2str(sigma,'%5.3f'),'.png'));
hold off;

plot((propx/2),propefield,'-r');
title(strcat('E_y versus x (in lambda), symmetric case, sigma=',num2str(sigma,'%5.3f')));
xlabel('x / $\lambda$', 'Interpreter', 'LaTex');
ylabel('E_y');
hold on;
plot((combx/2),combefield,'-b');
plot((darkx/2),darkefield,'-g');
legend('propagating only','combined', 'dark');
print('-dpng',strcat('plots/efieldoverlaylambdasigma',num2str(sigma,'%5.3f'),'.png'));
hold off;


end