%% reconstruct the LUT for speedy processing
function G0=get_G(isza,G)


solz=cell2mat({G(:).solz});
% solz=[0,15,30,45,60,75,80,88]; 

sidx=find(solz-isza>=0);
if isza==0
    isza1=solz(sidx(1));
    sid1=sidx(1);
else
    isza1=solz(sidx(1)-1);
    sid1=sidx(1)-1;
end

 G1=G(sid1).LUT;
 isza2=solz(sidx(1));
 sid2=sidx(1);
 G2=G(sid2).LUT;

if isza==0
    G0=G1;
else
    G0=(isza2-isza)/(isza2-isza1)*G1+(isza-isza1)/(isza2-isza1)*G2;
end


end
