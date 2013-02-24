function singdarkscript(secondinterface, thetamax, etacutoff)
gs='0.001';
data(:,:,1) = load(strcat('data/singlinedarkres',num2str(thetamax, '%3.1f'),'degs','3.1eta',gs,'sigmatilde',num2str(secondinterface,'%2.1f'),'secint',num2str(etacutoff,'%2i'),'etalimit.dat'));
x=data(:,1);
intensity=data(:,3);


plot(x,intensity,'-r');
title(strcat('Intensity profile (EE*) at image, thetamax: ', num2str(thetamax),' degs','kxcutoff',num2str(etacutoff,'%2i')));
xlabel('x ,d_s');
ylabel('Intensity (EE*)');
hold on;
%hilbertfit=abs(hilbert(ximageaxis));
%plot(xarray,hilbertfit,'-g');
%ylim([0 7e-03])
print('-dpng',strcat('plots/singlinedarkimageprofile','kxcutoff',num2str(etacutoff,'%2i'),num2str(thetamax),'degsthetamax.png'));
hold off;


%get fwhm
[maxval,maxindex]=max(intensity);
firsthalfx=x(1:maxindex);
secondhalfx=x(maxindex:length(x));
firsthalfint=intensity(1:maxindex);
secondhalfint=intensity(maxindex:length(intensity));

leftlowval = find(firsthalfint<=(maxval/2),1,'last');
lefthighval = find(firsthalfint>=(maxval/2),1,'first');
leftval=(firsthalfx(leftlowval)+firsthalfx(lefthighval))/2;

rightlowval = find(secondhalfint>=(maxval/2),1,'last');
righthighval = find(secondhalfint<=(maxval/2),1,'first');
rightval=(secondhalfx(rightlowval)+secondhalfx(righthighval))/2;

fwhm=rightval-leftval;

fid = fopen(strcat('data/singlinedarkplotdata',num2str(secondinterface),'secint.txt'), 'a');
outdata = [fwhm;thetamax;etacutoff; maxval];
fprintf(fid, '%6.4f %3.1f %3.1f %6.4e \n', outdata);
fclose(fid);



end