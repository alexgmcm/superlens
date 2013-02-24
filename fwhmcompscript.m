propdata(:,:,1) = load('data/singlinepropplotdata3secint.txt');
propfwhm=propdata(:,1);
logpropfwhm=log10(propfwhm);
propthetamax=propdata(:,2);
propintensity=propdata(:,4);

combdata(:,:,1)= load('data/singlinecombplotdata3secint.txt');
combfwhm=combdata(:,1);
logcombfwhm=log10(combfwhm);
combthetamax=combdata(:,2);
combintensity=combdata(:,4);

darkdata(:,:,1)=load('data/singlinedarkplotdata3secint.txt');
darkfwhm=darkdata(:,1);
logdarkfwhm=log10(darkfwhm);
darkthetamax=darkdata(:,2);
darkintensity=darkdata(:,4);

plot(propthetamax,propfwhm,'-r');
title('FWHM versus thetamax value (in degrees), symmetric case');
xlabel('thetamax / degrees');
ylabel('FWHM/source distance');
hold on;
plot(combthetamax,combfwhm,'-b');
plot(darkthetamax,darkfwhm,'-g');
legend('propagating only','combined', 'dark');
print('-dpng','plots/fwhmthetaplot.png');

hold off;


plot(propthetamax,propintensity,'-r');
title('Intensity versus thetamax value (in degrees), symmetric case');
xlabel('thetamax / degrees');
ylabel('Intensity (EE*)');
hold on;
plot(combthetamax,combintensity,'-b');
plot(darkthetamax,darkintensity,'-g');
legend('propagating only','combined', 'dark');
print('-dpng','plots/intensitythetaplot.png');
hold off

