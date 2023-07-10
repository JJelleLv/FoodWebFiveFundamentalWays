function [H]=func_HurDet(F)

%% nr of H to determine
NRH=length(F(:,1));

%% determine H
H=zeros(NRH,1);
Fvec=[zeros(NRH-2,1);1;-F;zeros(NRH+1,1)];
for HNR=1:NRH
    
    %% make matrix with feedbacks
    Hmat=zeros(HNR,HNR);
    for Hi=1:HNR
        Hmat(:,Hi)=flipud(Fvec([NRH-HNR+(Hi-1)*2+1:NRH+(Hi-1)*2],1));
    end
    
    %% Hurwitz determinant
    H(HNR,1)=det(Hmat);
    
end