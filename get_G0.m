%% relocate the G values for specific solar zenith angle
function G0=get_G0(isza,G_LUT)

solz=[0,15,30,45,60,75,80,88]; 
% [~,sidx]=min(abs(solz-isza));
% isza=solz(sidx);
sidx=find(solz-isza>=0);
if isza==0
    isza1=solz(sidx(1));
else
    isza1=solz(sidx(1)-1);
end
isza2=solz(sidx(1));

senz=[0,10,20,30,40,50,60,70,80,87.5];
phi=[0,15,30,45,60,75,90,105,120,135,150,165,180]; 
ntheta=length(senz);
nphi=length(phi);
G0=zeros(ntheta,nphi,4);
G1=G0;G2=G0;

for i = 1: ntheta
    for j=1:nphi
        if i==1
            G1(i,j,1:4)=G_LUT(1,4:7);     
            G2(i,j,1:4)=G_LUT(1,4:7);     
        else
                ivza=senz(i);
                ivaa=phi(j); 
                idx1=find(G_LUT(:,1)==isza1);
                idx2=find(G_LUT(:,2)==ivza);
                idx3=find(G_LUT(:,3)==ivaa);
                idx_temp=intersect(idx1,idx2);
                idx=intersect(idx_temp,idx3);
                G1(i,j,1:4)=squeeze(G_LUT(idx,4:7)); 
                
                idx11=find(G_LUT(:,1)==isza2);
                idx_temp2=intersect(idx11,idx2);
                idx_3=intersect(idx_temp2,idx3);
                G2(i,j,1:4)=squeeze(G_LUT(idx_3,4:7)); 
                
                
        end
    end
end

if isza==0
    G0=G1;
else
    G0=(isza2-isza)/(isza2-isza1)*G1+(isza-isza1)/(isza2-isza1)*G2;
end


end
