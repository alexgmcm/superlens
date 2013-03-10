function plotfwhmioscript(sigma)

cdata(:,:)=load(strcat('data/singlinecombinedplotdata',num2str(sigma,'%5.3f'),'3.0secint.txt'));
cfwhm=cdata(:,1);
thetamax=cdata(:,2);
cmaxval=cdata(:,3);

pdata(:,:)=load(strcat('data/singlineproponlyplotdata',num2str(sigma,'%5.3f'),'3.0secint.txt'));
pfwhm=pdata(:,1);
pthetamax=pdata(:,2);
pmaxval=pdata(:,3);

ddata(:,:)=load(strcat('data/singlinedarkonlyplotdata',num2str(sigma,'%5.3f'),'3.0secint.txt'));
dfwhm=ddata(:,1);
dmaxval=ddata(:,3);

plot(thetamax,cfwhm,'-r');
title(strcat('FWHM versus thetamax (in degrees), symmetric case, sigma: ', num2str(sigma,'%5.3f')) );
xlabel('thetamax / degrees');
ylabel('FWHM/source distance');
hold on;
plot(pthetamax,pfwhm,'-b');
plot(thetamax,dfwhm,'-g');
legend('combined','propagating only','dark');
print('-dpng',strcat('plots/sigma',num2str(sigma,'%5.3f'),'fwhmthetaplot.png'));

hold off;



plot(thetamax,cmaxval,'-r');
title(strcat('Intensity versus thetamax (in degrees), symmetric case, sigma: ', num2str(sigma,'%5.3f')) );
xlabel('thetamax / degrees');
ylabel('Intensity (EE*)');
hold on;
plot(pthetamax,pmaxval,'-b');
plot(thetamax,dmaxval,'-g');
legend('combined','propagating only','dark');
print('-dpng',strcat('plots/sigma',num2str(sigma,'%5.3f'),'intensitythetaplot.png'));

hold off;



end