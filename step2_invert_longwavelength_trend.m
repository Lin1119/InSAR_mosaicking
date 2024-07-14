
function [lt]=step2_invert_longwavelength_trend(G_all, D_all, in, s, q, meanll, lonlat_data)
% Step 2 - inverting long wavelength trend for individual InSAR track

% INPUTS:
% G_all ------ jumbo G matrix
% D_all ------ all data points
% in --------- parameter index
% q ---------- track index
% s ---------- standard deviation of velocities for all points
% lonlat_data --------- lon/lat coordinates of InSAR track

% OUTPUTS:
% lt ------------------ inverted long wavalength trend

% By Lin Shen -- University of Leeds

S=diag(s.^2);
m=lscov(G_all,D_all,S);
lonlat_data(:,1:2)=lonlat_data(:,1:2)-repmat(meanll,length(lonlat_data),1);
lt=m(in+q+1)*(lonlat_data(:,1))+m(in+q+2)*(lonlat_data(:,2))+m(in+q);




