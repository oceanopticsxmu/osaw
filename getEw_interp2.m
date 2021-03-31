function [Ew]=getEw_interp2(Lw,senz,phi) 

%Code by Xiaolong Yu (xlyu@xmu.edu.cn)

%% imput 
% Lw:       imput water-leaving radiance, matrix, with a dimension of ntheta*nphi
% senz:     the value of theta (viewing zenith angle) for each quad, vector (0-90), in degree
% dphij:    interval of the phi (viewing azimuth angle) in the two neighbouring quads, assuming evenly
%           distributed phi in the grid. dphij=360/nphi.
% nphi:     number of phi to partition the grid (vector, dimension, nphi*1), 
%           determines the angular resolution of the quad layout, delta_phi
% ntheta:   number of thetav to partition the grid (vector, dimension, ntheta*1), 
%           determines the angular resolution of the quad layout, delta_theta

%% output 
% Ew: upwelling irradiance contributed by water leaving signal
%     via integral of Lw over all directions
%     theta=[0 90], phi=[0 360].

%% main 

% [ntheta,nphi]=size(Lw);    
% dphij=max(phi)/(nphi-1); 

% if size(thetav,2)~=1
%     thetav=thetav';
% end
% imput data quality check 
% if ntheta ~= length(senz)
%     disp('Warning: radiance and viewing angles do not match!');
%     error('The dimension of input radiance does not match the input viewing angles');
% end
  
deg2rad=180/pi();        % deg to radian
senz=senz/deg2rad;   % convert to radian
phi=phi/deg2rad;

Lw_temp=Lw;               % w_temp is the temporary index for the computed irradiance at each quad 
Lw_temp(Lw_temp<0.00001)=0.00001;   % avoid negative radiance inputs 

% intergal for an interval of 1 deg

step=1/deg2rad;
% senz_intp=senz(1):step:senz(end);  % range: 0~ 87.5
% senz_intp(end)=senz(end);

senz_intp=0:step:pi()/2;    % extend to 90 (pi/2) 
phi_intp=phi(1):step:phi(end);
phi_intp(end)=phi(end);

[X1,Y1]=meshgrid(phi,senz);
[X2,Y2] = meshgrid(phi_intp,senz_intp);  % integral

Lw_intp = interp2(X1,Y1,Lw_temp,X2,Y2,'spline');

cos_temp=cos(senz_intp');
sin_temp=sin(senz_intp'); 

itemp=Lw_intp.*cos_temp.*sin_temp;

Ew0_temp=trapz(senz_intp,itemp,1);  % first integral over theta
Ew=trapz(phi_intp,Ew0_temp);

% Ew0_temp=trapz(phi,itemp,2);  % first integral over phi
% Ew=trapz(thetav,Ew0_temp);
Ew=2*Ew;  % for phi from pi to 2pi

end