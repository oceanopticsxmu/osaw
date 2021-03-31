function [a,bbp] = qaa_Lee11(wl,Rrs,G,aw,bbw)
% Combination of QAA_v6 and that IOPs inversion algroithm proposed in Lee
% et al., 2011
% coded by Xiaolong Yu (xlyu@xmu.edu.cn), Feb 20,2021

% Inputs:
%   Rrs: remote-sensing reflectance
%   wl: wavelength
%   G: G values for specific angular geometry, corresping to Rrs angle info
%   aw: pure seawater absorption coefficients 
%   bbw: pure seawater backscattering coefficients

% Outputs:
%   a, bb: total absorbtion and backscattering
%   anw: nonwater absorption, anw= a - aw
%   bbp: particulate backscattering
%   aph: phytoplankton abs
%   acdm: CDM absorption

% *************************************************************
% Citation: Lee Z. P., IOCCG, online link, quasi-analytical algorithm (QAA-v6); not published
% Lee et al., 2011. An inherent-optical-property-centered approach to correct the angular effects in water-leaving radiance
% ***************************************************************************
 
nwl=length(wl);
rrs(1:nwl) = Rrs./(0.52+1.7*Rrs);
% % id412 = find(abs(wl-412)==min(abs(wl-412))); 
id443 = find(abs(wl-443)==min(abs(wl-443))); 
id490 = find(abs(wl-490)==min(abs(wl-490))); 
id555 = find(abs(wl-555)==min(abs(wl-555))); 
id670 = find(abs(wl-670)==min(abs(wl-670)));  
id555=id555(1);
id670=id670(1);


G0w=G(1);
G1w=G(2); 
G0p=G(3);
G1p=G(4);

if Rrs(id670) < 0.0015     
    
    wl_ref = wl(id555);
    id_ref = id555;     
    %%%%
    ki = log10((rrs(id443)+rrs(id490))/(rrs(id_ref) + 5*rrs(id670)*rrs(id670)/rrs(id490)));
    a_ref = aw(id_ref) + 10.^(-1.146 - 1.366 * ki - 0.469*ki.^2);
else
    wl_ref = 670;
    id_ref = id670;
    a_ref = aw(id_ref) + 0.39 * (rrs(id670)/(rrs(id443) + rrs(id490)))^1.14;
end


bbw_ref=bbw(id_ref);
Rrs_ref=Rrs(id_ref);

A=G0p+G1p-Rrs_ref;
B=G0w*bbw_ref+G0p*(a_ref+bbw_ref)-2*(a_ref+bbw_ref)*Rrs_ref;
C=G0w*(a_ref+bbw_ref)*bbw_ref-Rrs_ref*(a_ref+bbw_ref).^2+G1w*bbw_ref.^2;

bbp_ref=(sqrt(B.^2-4*A*C)-B)/(2*A);

eta= 2.0*(1-1.2*exp(-0.9*rrs(id443)/rrs(id555))); 
bbp(1:nwl)=bbp_ref *(wl_ref./wl).^eta;

% for i=1:nwl
    X=G0w*bbw+G0p*bbp;
    Y=G1w*bbw.^2+G1p*bbp.^2;
    kappa=(sqrt(X.^2+4*Rrs.*Y)+X)./(2*Rrs);
% end

a(1:nwl) =kappa-bbp-bbw;


end


