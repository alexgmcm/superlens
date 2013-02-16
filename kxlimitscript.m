data(:,:,1) = load('data/oldsinglinecombplotdata3secint.txt');
fwhm=data(:,1);
etacutoff=data(:,3);

plot(etacutoff,fwhm,'-r');
title('FWHM versus kx cut-off value (in terms of eta), symmetric case');
xlabel('kx cut-off , eta');
ylabel('FWHM, source distance');
hold on;
%hilbertfit=abs(hilbert(ximageaxis));
%plot(xarray,hilbertfit,'-g');
%ylim([0 7e-03])
print('-dpng','plots/symmetrickxcutoffplot.png');
hold off;