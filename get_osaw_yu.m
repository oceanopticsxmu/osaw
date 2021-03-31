function [osaw,a,bb]=get_osaw_yu(wl,Rrs,isza,opt)
% Code by Xiaolong Yu (xlyu@xmu.edu.cn) @ 2021-02-03
% revision history: 
% reference,Yu et al., 2021, Estimation water-leaving albedo from ocean color remote sensing
%% imput
% wl:            wavelength (nm)
% Rrs:           (satelite)measured remote sensing reflectance, in sr-1 
% opt:          structure for the input, 
% opt.G_LUT:    Look-up-table for G values describe the angular
%               distribution of Rrs (Lee et al., 2011)
% opt.inverse:  inverse algorithms for IOPs, default 'QAA'
% opt.aw: pure seawater absortion, either 'Lee'  (Lee et al., 2015) or
%             default (Pope & Fry, 1997), optional 
% opt.bbw: pure seawater backscattering, either 'Zhang'  (Zhang et al., 2009) or
%             default (Morel, 1974), optional 

%% output 
% osaw: water leaving ocean surface albedo at a specific wavelength at wl 
% functions to call:

%% main code 
n_wl=length(wl);
wl(1:n_wl)=wl;
Rrs(1:n_wl)=Rrs; 

% angular quad-->  Hydrolight simulation
senz=[0,10,20,30,40,50,60,70,80,87.5];
phi=[0,15,30,45,60,75,90,105,120,135,150,165,180]; 
ntheta=length(senz);
nphi=length(phi);

if exist('opt') == 0 
    opt = struct; 
    'a structure of input is required';
    return
end  
%% select aw and bbw 
    aw(1:n_wl) =opt.aw;
    bbw(1:n_wl) = opt.bbw;    

%% IOPs inverse algorithm
    G=opt.G_ref;
if isfield(opt,'inverse') == 0  
    'an inverse algroithm is required';
    return
elseif strcmp(opt.inverse,'QAA') 
    [a,bbp] = qaa_Lee11(wl,Rrs,G,aw,bbw);  % use QAA derived IOPs 
else
    [a,bbp] = qaa_Lee11(wl,Rrs,G,aw,bbw);  % set as default, may be replaced later
  
end

if isfield(opt,'bbp') == 1  
   bbp(1:n_wl)=opt.bbp;
end

if isfield(opt,'a') == 1  
   a(1:n_wl)=opt.a;
end

bb=bbp+bbw;
k=a+bb; % kappa
%% get G values and compute angular Rrs
G0=get_G(isza,opt.G_LUT); 

Angular_Rrs=zeros(ntheta,nphi,n_wl);
  
for i = 1: ntheta
    for j=1:nphi
       Angular_Rrs(i,j,1:n_wl)=(G0(i,j,1)+G0(i,j,2)*bbw./k).*bbw./k+(G0(i,j,3)+G0(i,j,4)*bbp./k).*bbp./k; 
    end
end

for k = 1:length(wl)
    A_Rrs=Angular_Rrs(:,:,k);
    osaw(k)=getEw_interp2(A_Rrs,senz,phi);
end


end
