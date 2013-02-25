propdata(:,:,1) = load('data/proponlyimageline90.0degs3.1eta0.001sigmatilde3.0secint30etalimit.dat');
propx=propdata(:,1);
propintensity=propdata(:,3);
propefield=propdata(:,4);

combdata(:,:,1) = load('data/combinedimageline90.0degs3.1eta0.001sigmatilde3.0secint30etalimit.dat');
combx=combdata(:,1);
combintensity=combdata(:,3);
combefield=combdata(:,4);

darkdata(:,:,1) = load('data/darkonlyimageline90.0degs3.1eta0.001sigmatilde3.0secint30etalimit.dat');
darkx=darkdata(:,1);
darkintensity=darkdata(:,3);
darkefield=darkdata(:,4);

plot(propx,propintensity,'-r');
title('Intensity versus x (in dsource), symmetric case, sigma=0.001');
xlabel('x / dsource');
ylabel('Intensity , EE*');
hold on;
plot(combx,combintensity,'-b');
plot(darkx,darkintensity,'--g');
legend('propagating only','combined', 'dark');
print('-dpng','plots/intensityoverlaysigma0.001.png');
hold off;

plot(propx,propefield,'-r');
title('E_y versus x (in dsource), symmetric case, sigma=0.001');
xlabel('x / dsource');
ylabel('E_y');
hold on;
plot(combx,combefield,'-b');
plot(darkx,darkefield,'--g');
legend('propagating only','combined', 'dark');
print('-dpng','plots/efieldoverlaysigma0.001.png');
hold off;