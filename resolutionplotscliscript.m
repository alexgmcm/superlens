function resolutionplotscliscript(thetamax)

indata= load(strcat('data/plotdata',num2str(thetamax),'degs.txt'));
slabthicknessvals=indata(:,1);
imagepositions=indata(:,2);
fwhmvals=indata(:,3);
expectedpendryvals=indata(:,4);


plot(expectedpendryvals,imagepositions,'-xr')
title(strcat('Test of Pendry Lens Equation, thetamax: ', num2str(thetamax),' degs'));
xlabel('Expected image distance from source');
ylabel('Actual image distance from source');
print('-dpng',strcat('plots/pendrytest',num2str(thetamax),'degsthetamax.png'));



plot(slabthicknessvals,fwhmvals,'-xr')
title(strcat('Resolution for given slab thickness, thetamax: ', num2str(thetamax),' degs'));
xlabel('Slab thickness');
ylabel('FWHM of intensity of image');
print('-dpng',strcat('plots/resplot',num2str(thetamax),'degsthetamax.png'));


end