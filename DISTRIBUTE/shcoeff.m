function [sh_coeff, diff] = shcoeff(points,ADC_allcmpts_alldir)

% plot ngdir original directions in which the HADC was simulated and 900 interpolated directions
% 
% Input: 
%     1. points (ngdir directions)
%     2. ADC_allcmpts_alldir
%     
% Output:
%     1 figure with title of "ADC in ngdir directions"

YY = spherical_harmonics(points(:,1),points(:,2),points(:,3),ones(size(points,1),1));
sh_coeff = YY\ADC_allcmpts_alldir;

[sph_pts, ~] = spheresurface_regularpoints(1,900);   
YY_sph = spherical_harmonics(sph_pts(:,1),sph_pts(:,2),sph_pts(:,3),ones(size(sph_pts,1),1));
ADC_interp = YY_sph*sh_coeff;
ADC_alldir = ADC_interp.*sph_pts;

location = abs(points*(sph_pts')-1)<3e-3;
len = length(points);
ind = ones(len, 1);
for i= 1:len
    a = find(location(i, :));
    ind(i) = a(1);
end
diff = ADC_allcmpts_alldir-ADC_interp(ind);
