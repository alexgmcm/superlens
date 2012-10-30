function [ refractedAngle ] = snells( incidentAngle )
%snells Returns expected angle of refraction from snells law for eps2=10
incidentAngle = incidentAngle *(pi/180);
refractedAngle=asin(sin(incidentAngle)./sqrt(10))*(180/pi);


end

