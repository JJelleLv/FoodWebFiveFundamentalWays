function [INT_list_c,NRint_c]=func_relInt_c(NRspec,NRBASAL,TL_list,MIN_C,MAX_C)

% make competitive interaction list, speciesNR, speciesNR, strength_k
IntNR_c=0;
INT_list_c=[];
if NRBASAL>1
    INT_list_c=nan((nchoosek(NRBASAL,2)+NRBASAL),3);
    for SpecNR_i=1:NRspec
        if TL_list(SpecNR_i,1)==1;
            
            inter_c=1;
            
            IntNR_c=IntNR_c+1;
            INT_list_c(IntNR_c,[1:3])=[SpecNR_i SpecNR_i inter_c];
            
            for SpecNR_j=SpecNR_i+1:NRspec
                if TL_list(SpecNR_j,1)==1;
                    
                    intra_c=rand.*(MAX_C-MIN_C)+MIN_C; %% could be a random number
                    
                    IntNR_c=IntNR_c+1;
                    INT_list_c(IntNR_c,[1:3])=[SpecNR_i SpecNR_j intra_c];
                    
                end
            end
        end
    end
end
NRint_c=IntNR_c;