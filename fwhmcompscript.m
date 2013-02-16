propdata(:,:,1) = load('data/plotdata3secint.txt');
propfwhm=propdata(:,3);
logpropfwhm=log10(propfwhm);
propthetamax=propdata(:,5);

combdata(:,:,1)= load('data/singlinecombplotdata3secint.txt');
combfwhm=combdata(:,1);
logcombfwhm=log10(combfwhm);
combthetamax=combdata(:,2);

plot(propthetamax,propfwhm,'-r');
title('log base 10 FWHM versus thetamax value (in degrees), symmetric case');
xlabel('thetamax , degrees');
ylabel('log base 10 of FWHM, source distance');
hold on;
plot(combthetamax,combfwhm,'-b');
%hilbertfit=abs(hilbert(ximageaxis));
%plot(xarray,hilbertfit,'-g');
%ylim([0 7e-03])
print('-dpng','plots/log10fwhmthetaplot.png');
legend('propagating only','combined');
hold off;