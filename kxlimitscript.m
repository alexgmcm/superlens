data(:,:,1) = load('data/oldsinglinecombplotdata3secint.txt');
fwhm=data(:,1);
etacutoff=data(:,3);

plot(etacutoff,(fwhm/2),'-r','LineWidth',4);
title('FWHM versus kx cut-off value (in terms of eta), symmetric case','FontSize',16);
xlabel('kx cut-off / eta','FontSize',16);
ylabel('FWHM / lambda','FontSize',16);
hold on;
%hilbertfit=abs(hilbert(ximageaxis));
%plot(xarray,hilbertfit,'-g');
%ylim([0 7e-03])
print('-depsc2','plots/symmetrickxcutoffplot.eps');
hold off;