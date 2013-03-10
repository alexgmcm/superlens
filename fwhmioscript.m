function fwhmioscript(thetamax,sigma,string)

data(:,:,1) = load(strcat('data/',string,'imageline',num2str(thetamax,'%4.1f'),'degs3.1eta',num2str(sigma,'%5.3f'),'sigmatilde3.0secint30etalimit.dat'));
x=data(:,1);
intensity=data(:,3);
efield=data(:,4);

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

fid = fopen(strcat('data/singline',string,'plotdata',num2str(sigma,'%5.3f'),'3.0secint.txt'), 'a');
outdata = [fwhm;thetamax; maxval];
fprintf(fid, '%6.4f %3.1f %6.4e \n', outdata);
fclose(fid);












end