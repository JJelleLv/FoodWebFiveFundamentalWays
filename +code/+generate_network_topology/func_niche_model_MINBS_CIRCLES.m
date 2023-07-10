function [NETnichevalue_uns,NETnichemax_uns,NETnichemin_uns,NETcvalue_uns,NETrvalue_uns,NET_uns, ...
    NRBASAL,MINNRBASAL,INT_list_NET_uns,BASAL_list,TOPPRED_list,TL_list,NRint_uns,bniche, ...
    CIRCLE_Int_list,CIRCLE_Int_list_NR,CIRCLENR,CIRCLE_nodes,CIRCLE_links,Spec_NRprey_count] = func_niche_model_MINBS_CIRCLES(S,CONN,dCONN,MINNRBASAL)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Function that generates food web topology %%%%%%%%%%%%%%%%%
%%%% according to niche model (Williams and Martinez, 2000) %%%%
%%%% Code by Jelle Lever, 2014 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% updated by Adrian Etten %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Generate networks %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while true
    while true
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Generate network with niche model %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % initialisation - build foodNET
        NETnichevalue_uns=sort(rand(S,1));
        bniche=1/(2*CONN)-1;
        NETrvalue_uns=(betarnd(1,bniche,S,1)).*NETnichevalue_uns;
        NETcvalue_uns=(NETnichevalue_uns-(NETrvalue_uns/2)).*rand(S,1)+(NETrvalue_uns/2);
        NETnichemin_uns=NETcvalue_uns-(NETrvalue_uns/2);
        NETnichemax_uns=NETcvalue_uns+(NETrvalue_uns/2);
        
        NET_uns=zeros(S,S);
        for PREY=1:S
            for PRED=1:S
                if NETnichevalue_uns(PREY,1) > NETnichemin_uns(PRED,1) && NETnichevalue_uns(PREY,1) < NETnichemax_uns(PRED,1)
                    NET_uns(PREY,PRED)=1;
                end
            end
        end
        
        %% do not allow when the not the desired NR of basal species
        BASAL_list=(sum(NET_uns)'==0);
        NRBASAL=sum(BASAL_list);
        TOPPRED_list=(sum(NET_uns,2)==0);
        
        %% do not allow basal species without any links
        TOPBASALnr=sum(BASAL_list.*TOPPRED_list);
        
        %% do not allow when dCONN outside range
        real_dCONN=abs(CONN-(sum(sum(NET_uns))/(S^2)));
        
        %if NRBASAL_outp==NRBASAL && TOPBASALnr==0
        if NRBASAL>=MINNRBASAL && TOPBASALnr==0 && real_dCONN<=dCONN
            break
        end
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% CHECK TROPHIC SIMILARITY %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    TROPHSIM=0;
    for Spec_i=1:S
        for Spec_j=(Spec_i+1):S
            if (NET_uns(:,Spec_i)==NET_uns(:,Spec_j)) %% the same PREY
                if (NET_uns(Spec_i,:)==NET_uns(Spec_j,:)) %% the same PREDATOR
                    %Spec_i
                    %Spec_j
                    TROPHSIM=TROPHSIM+1;
                    %stop
                end
            end
        end
    end
    
    NOOVERLAP=0;
    %     for Spec_i=1:S
    %         for Spec_j=(Spec_i+1):S
    %             if (NET_uns(:,Spec_i)==NET_uns(:,Spec_j)) %% the same PREY
    %                 if BASAL_list(Spec_i,1)==0 && BASAL_list(Spec_j,1)==0 %% not BASAL
    %                     %Spec_i
    %                     %Spec_j
    %                     NOOVERLAP=1;
    %                     %stop
    %                 end
    %             end
    %         end
    %     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% DETERMINE CIRLES, TROPHIC LEVEL and so on ... %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% only the do the next steps when no similarity
    if TROPHSIM==0 && NOOVERLAP==0
        
        %% make interaction list
        IntNR_uns=0;
        INT_list_NET_uns=[];
        for PREY=1:S
            for PRED=1:S
                if NET_uns(PREY,PRED)==1;
                    IntNR_uns=IntNR_uns+1;
                    INT_list_NET_uns(IntNR_uns,1)=PREY;
                    INT_list_NET_uns(IntNR_uns,2)=PRED;
                end
            end
        end
        NRint_uns=IntNR_uns;
        
        %% find Circles in the network (not loops)
        CIRCLE_nodes={};
        CIRCLE_links={};
        
        CIRCLENR=0;
        CIRCLE_Int_list=[];
        CIRCLE_Int_list_NR=0;
        
        for SpecNR=1:S
            
            DIRECTED_BRANCH_nodes={};
            DIRECTED_BRANCH_links={};
            DIRECTED_BRANCH_finalnodes={};
            
            if (sum(INT_list_NET_uns(:,2)==SpecNR))>=1
                
                %% starting point branches
                BRANCHNR_per_DIRECTED_BRANCH_length=zeros(S,1);
                BRANCHNR_per_DIRECTED_BRANCH_length(1,1)=1;
                DIRECTED_BRANCH_nodes{1,1}=SpecNR;
                DIRECTED_BRANCH_links{1,1}=[];
                DIRECTED_BRANCH_finalnodes{1,1}=SpecNR;
                
                %% go through all branches with length... till S
                for DIRECTED_BRANCH_length=2:S
                    
                    for NRBRANCH=1:BRANCHNR_per_DIRECTED_BRANCH_length((DIRECTED_BRANCH_length-1),1)
                        
                        %% go through interaction list
                        
                        for IntNR_t=1:NRint_uns
                            
                            PREY=INT_list_NET_uns(IntNR_t,1);
                            PRED=INT_list_NET_uns(IntNR_t,2);
                            
                            if PRED~=PREY
                                
                                NODEONBRANCHPREY=sum(DIRECTED_BRANCH_nodes{NRBRANCH,(DIRECTED_BRANCH_length-1)}==PREY);
                                if DIRECTED_BRANCH_finalnodes{NRBRANCH,(DIRECTED_BRANCH_length-1)}==PRED && NODEONBRANCHPREY==0
                                    
                                    BRANCHNR_per_DIRECTED_BRANCH_length(DIRECTED_BRANCH_length,1)=BRANCHNR_per_DIRECTED_BRANCH_length(DIRECTED_BRANCH_length,1)+1;
                                    DIRECTED_BRANCH_finalnodes{BRANCHNR_per_DIRECTED_BRANCH_length(DIRECTED_BRANCH_length,1),DIRECTED_BRANCH_length}=PREY;
                                    
                                    DIRECTED_BRANCH_nodes{BRANCHNR_per_DIRECTED_BRANCH_length(DIRECTED_BRANCH_length,1),DIRECTED_BRANCH_length}=DIRECTED_BRANCH_nodes{NRBRANCH,(DIRECTED_BRANCH_length-1)};
                                    DIRECTED_BRANCH_nodes{BRANCHNR_per_DIRECTED_BRANCH_length(DIRECTED_BRANCH_length,1),DIRECTED_BRANCH_length}(1,DIRECTED_BRANCH_length)=PREY;
                                    
                                    DIRECTED_BRANCH_links{BRANCHNR_per_DIRECTED_BRANCH_length(DIRECTED_BRANCH_length,1),DIRECTED_BRANCH_length}=DIRECTED_BRANCH_links{NRBRANCH,(DIRECTED_BRANCH_length-1)};
                                    DIRECTED_BRANCH_links{BRANCHNR_per_DIRECTED_BRANCH_length(DIRECTED_BRANCH_length,1),DIRECTED_BRANCH_length}(1,(DIRECTED_BRANCH_length-1))=IntNR_t;
                                    
                                end
                            end
                            
                            %% Circle found!
                            if PREY==SpecNR && DIRECTED_BRANCH_finalnodes{NRBRANCH,(DIRECTED_BRANCH_length-1)}==PRED
                                
                                CIRCLENR=CIRCLENR+1;
                                
                                CIRCLE_nodes{CIRCLENR,1}=DIRECTED_BRANCH_nodes{NRBRANCH,(DIRECTED_BRANCH_length-1)};
                                CIRCLE_nodes{CIRCLENR,1}(1,DIRECTED_BRANCH_length)=PREY;
                                
                                CIRCLE_links{CIRCLENR,1}=DIRECTED_BRANCH_links{NRBRANCH,(DIRECTED_BRANCH_length-1)};
                                CIRCLE_links{CIRCLENR,1}(1,(DIRECTED_BRANCH_length-1))=IntNR_t;
                                
                                %% go through links on Circle and add them to C
                                for NR_Int_CIRCLE=1:(DIRECTED_BRANCH_length-1)
                                    Int_CIRCLE=CIRCLE_links{CIRCLENR,1}(1,NR_Int_CIRCLE);
                                    %% if not on lis
                                    if (sum(CIRCLE_Int_list==Int_CIRCLE))==0
                                        CIRCLE_Int_list_NR=CIRCLE_Int_list_NR+1;
                                        CIRCLE_Int_list(CIRCLE_Int_list_NR,1)=Int_CIRCLE;
                                    end
                                end
                                
                            end
                            
                        end
                    end
                end
            end
        end
        
        CIRCLE_Int_list=sort(CIRCLE_Int_list);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Calculate trophic level - connected to basal %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        TL_list=double(BASAL_list);
        FOUND_list=double(BASAL_list);
        
        %% find number pred and nr prey per spec
        Spec_NRprey_count=zeros(S,1);
        for SpecNR=1:S
            Spec_NRprey_count(SpecNR,1)=sum(INT_list_NET_uns(:,2)==SpecNR);
        end
        
        CHECKFOUND=0;
        while CHECKFOUND==0
            
            Spec_NRprey_count_nocircle=Spec_NRprey_count;
            
            COUNT_prey_found=zeros(S,1);
            SUM_TL_PREY=zeros(S,1);
            for IntNR_uns=1:NRint_uns
                
                PREY=INT_list_NET_uns(IntNR_uns,1);
                PRED=INT_list_NET_uns(IntNR_uns,2);
                
                if (sum(CIRCLE_Int_list==IntNR_uns))==0  %%PREY~=PRED %% shoulde become: interaction not on circle
                    
                    if FOUND_list(PREY,1)==1
                        COUNT_prey_found(PRED,1)=COUNT_prey_found(PRED,1)+1;
                        SUM_TL_PREY(PRED,1)=SUM_TL_PREY(PRED,1)+TL_list(PREY,1);
                    end
                    
                elseif (sum(CIRCLE_Int_list==IntNR_uns))>0 %%PREY==PRED %% shoulde become: interaction on circle
                    
                    Spec_NRprey_count_nocircle(PRED,1)=Spec_NRprey_count_nocircle(PRED,1)-1;
                    
                end
                
            end
            
            NEWFOUND=0;
            for SpecNR=1:S
                if COUNT_prey_found(SpecNR,1)==Spec_NRprey_count_nocircle(SpecNR,1) && FOUND_list(SpecNR,1)==0
                    
                    %(SUM_TL_PREY(SpecNR,1)/Spec_NRprey_count_noself(SpecNR,1))+1
                    TL_list(SpecNR,1)=(SUM_TL_PREY(SpecNR,1)/Spec_NRprey_count_nocircle(SpecNR,1))+1;
                    FOUND_list(SpecNR,1)=1;
                    NEWFOUND=1;
                end
            end
            
            if NEWFOUND==0 %% some species are disconnected from basal
                CHECKFOUND=-9999;
            elseif NEWFOUND==1
                SUMFOUND_list=sum(FOUND_list);
            end
            
            if SUMFOUND_list==S
                CHECKFOUND=1;
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%
        %%%%% Save data %%%%%
        %%%%%%%%%%%%%%%%%%%%%
        
        %% check if all connected to BASAL
        if (sum((BASAL_list+TL_list)>1))==S && CHECKFOUND==1 %% && TROPHSIM==0 && NOOVERLAP==0
 
            %% display status
            %fprintf('Niche Model generated')
            
            %% network found, so return
            return
            
        end
    end
end