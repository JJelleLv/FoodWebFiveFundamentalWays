function [F,Fpos,Fneg]=func_feedback_posneg(Jac)

%% NRspec
NRspec=length(Jac(1,:));

%% calculate feedback
F=zeros(NRspec,1);
Fpos=zeros(NRspec,1);
Fneg=zeros(NRspec,1);
for lLength=1:NRspec
    
    %% get all species combinations (stored to speedup large networks)    
    %allCombData=load(sprintf('results%sALL_comb_S%d%sALL_comb_%d',filesep,NRspec,filesep,lLength));
    %ALL_comb=allCombData.ALL_comb;
    %NRcomb=allCombData.NRcomb;
    ALL_comb=nchoosek([1:NRspec],lLength);
    NRcomb=length(ALL_comb(:,1));
    
    for COMBNR=1:NRcomb
        
        %fprintf('%d - %d/%d\n', lLength, COMBNR, NRcomb)
        
        %% determinant of each species combination
        COMB=ALL_comb(COMBNR,:);
        Jac_COMB=Jac(COMB,COMB);
        
        if mod(lLength,2)==1 %% change sign when odd
            F_COMB=det(Jac_COMB);
        else
            F_COMB=-det(Jac_COMB);
        end

        %% sum off feedback strengths
        F(lLength,1)=F(lLength,1)+F_COMB;
        if F_COMB>=0
            Fpos(lLength,1)=Fpos(lLength,1)+F_COMB;
        else
            Fneg(lLength,1)=Fneg(lLength,1)+F_COMB;
        end
    end

end