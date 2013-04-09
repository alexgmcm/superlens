function plotfwhmioscript(sigma)

cdata(:,:)=load(strcat('data/singlinecombinedplotdata',num2str(sigma,'%5.3f'),'3.0secint.txt'));
cfwhm=cdata(:,1);
thetamax=cdata(:,2);
width=thetatowidth(thetamax);
cmaxval=cdata(:,3);

pdata(:,:)=load(strcat('data/singlineproponlyplotdata',num2str(sigma,'%5.3f'),'3.0secint.txt'));
pfwhm=pdata(:,1);
%pthetamax=pdata(:,2);
pmaxval=pdata(:,3);

ddata(:,:)=load(strcat('data/singlinedarkonlyplotdata',num2str(sigma,'%5.3f'),'3.0secint.txt'));
dfwhm=ddata(:,1);
dmaxval=ddata(:,3);

plot(thetamax,cfwhm,'-r');
title(strcat('FWHM versus thetamax (in degrees), symmetric case, sigma: ', num2str(sigma,'%5.3f')) );
xlabel('thetamax / degrees');
ylabel('FWHM/source distance');
hold on;
plot(thetamax,pfwhm,'-b');
plot(thetamax,dfwhm,'-g');
legend('combined','propagating only','dark');
print('-dpng',strcat('plots/sigma',num2str(sigma,'%5.3f'),'fwhmthetaplot.png'));

hold off;


plot(width,(cfwhm/2),'-r','LineWidth',4);
title(strcat('FWHM versus width (in lambda), symmetric case, sigma: ', num2str(sigma,'%5.3f')) ,'FontSize',16);
xlabel(' \textrm{width} / $\lambda$','Interpreter','LaTex','FontSize',16);
ylabel('FWHM/ lambda','FontSize',16);
hold on;
plot(width,(pfwhm/2),'-b','LineWidth',4);
plot(width,(dfwhm/2),'--g','LineWidth',4);
leg=legend('combined','propagating only','dark');
set(leg,'FontSize',16);
xlim([0 30]);
print('-depsc2',strcat('plots/sigma',num2str(sigma,'%5.3f'),'fwhmwidthplot.eps'));

hold off;


plot(thetamax,cmaxval,'-r');
title(strcat('Intensity versus thetamax (in degrees), symmetric case, sigma: ', num2str(sigma,'%5.3f')) );
xlabel('thetamax / degrees');
ylabel('Intensity (EE*)');
hold on;
plot(thetamax,pmaxval,'-b');
plot(thetamax,dmaxval,'-g');
legend('combined','propagating only','dark');
print('-dpng',strcat('plots/sigma',num2str(sigma,'%5.3f'),'intensitythetaplot.png'));

hold off;

plot(width,cmaxval,'-r','LineWidth',4);
title(strcat('Intensity versus width (in lambda), symmetric case, sigma: ', num2str(sigma,'%5.3f')),'FontSize',16 );
xlabel(' \textrm{width} / $\lambda$','Interpreter','LaTex','FontSize',16);
ylabel('Intensity (EE*)','FontSize',16);
hold on;
plot(width,pmaxval,'-b','LineWidth',4);
plot(width,dmaxval,'-g','LineWidth',4);
leg=legend('combined','propagating only','dark');
set(leg,'FontSize',16);
xlim([0 30]);
print('-depsc2',strcat('plots/sigma',num2str(sigma,'%5.3f'),'intensitywidthplot.eps'));

hold off;



end