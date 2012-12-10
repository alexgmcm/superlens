function resolutionplotscliscript(secondinterface)

indata= load(strcat('data/plotdata',num2str(secondinterface),'secint.txt'));
slabthicknessvals=indata(:,1);
imagepositions=indata(:,2);
fwhmvals=indata(:,3);
expectedpendryvals=indata(:,4);
thetamaxvals=indata(:,5);


% plot(expectedpendryvals,imagepositions,'-xr');
% title(strcat('Test of Pendry Lens Equation, thetamax: ', num2str(thetamax),' degs'));
% xlabel('Expected image distance from source');
% ylabel('Actual image distance from source');
% print('-dpng',strcat('plots/pendrytest',num2str(thetamax),'degsthetamax.png'));



plot(thetamaxvals,fwhmvals,'-xr');
title(strcat('Resolution for given thetamax, second interface: ', num2str(secondinterface),' d_s'));
xlabel('thetamax, degrees');
ylabel('FWHM of intensity of image, d_s');
print('-dpng',strcat('plots/resplot',num2str(secondinterface),'secint.png'));


end