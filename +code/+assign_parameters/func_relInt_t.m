function [INT_list_t_REL]=func_relInt_t(NRspec,NRint_t,NET_uns,INT_list_NET_uns,NETnichevalue_uns,NETcvalue_uns,TL_list,Spec_NRprey_count,FRAC_int_nrnd,VAR_int,FRAC_fe_herb,MIN_fe_herb,MAX_fe_herb,MEAN_fe_herb,VAR_fe_herb,FRAC_fe_carn,MIN_fe_carn,MAX_fe_carn,MEAN_fe_carn,VAR_fe_carn,MIN_delta_herb,MAX_delta_herb,MIN_delta_carn,MAX_delta_carn,STRENGTH_SORTING)

%% take relative interaction strengths from dirichlet distribution
A_gamma=-(-1+4.*VAR_int+2.*FRAC_int_nrnd-(FRAC_int_nrnd.^2))./(8.*VAR_int);

Strth_list=cell(1,NRspec);
PREY_list=cell(1,NRspec);
for SpecNR=1:NRspec
    
    DEGREE=Spec_NRprey_count(SpecNR,1);
    if DEGREE>=1
        
        MIN_int=FRAC_int_nrnd./DEGREE;
        a_gamma_MAT=ones(1,DEGREE).*A_gamma;
        
        DIST_DIRI=gamrnd(a_gamma_MAT,1,1,DEGREE);
        DIST_DIRI=DIST_DIRI./sum(DIST_DIRI);
        DIST_DIRI=MIN_int+DIST_DIRI.*((1-FRAC_int_nrnd));
        
        Strth_list{SpecNR}=DIST_DIRI;
        
    end
    
    %% added for niche distribution of int strengths
    PREY_list{SpecNR}=find(NET_uns(:,SpecNR));
    
end

%% generate delta's and fe's
%% calculate a and b of beta dist
a_beta_fe_herb=((-MEAN_fe_herb.*VAR_fe_herb)+MIN_fe_herb.*VAR_fe_herb-MEAN_fe_herb.^3+2.*(MEAN_fe_herb.^2).*MIN_fe_herb-(MIN_fe_herb.^2).*MEAN_fe_herb+MAX_fe_herb.*(MEAN_fe_herb.^2)-2.*MAX_fe_herb.*MIN_fe_herb.*MEAN_fe_herb+MAX_fe_herb.*(MIN_fe_herb.^2))./((MAX_fe_herb-MIN_fe_herb).*VAR_fe_herb);
b_beta_fe_herb=(a_beta_fe_herb.*(-MEAN_fe_herb+MAX_fe_herb)/(MEAN_fe_herb-MIN_fe_herb));

a_beta_fe_carn=((-MEAN_fe_carn.*VAR_fe_carn)+MIN_fe_carn.*VAR_fe_carn-MEAN_fe_carn.^3+2.*(MEAN_fe_carn.^2).*MIN_fe_carn-(MIN_fe_carn.^2).*MEAN_fe_carn+MAX_fe_carn.*(MEAN_fe_carn.^2)-2.*MAX_fe_carn.*MIN_fe_carn.*MEAN_fe_carn+MAX_fe_carn.*(MIN_fe_carn.^2))./((MAX_fe_carn-MIN_fe_carn).*VAR_fe_carn);
b_beta_fe_carn=(a_beta_fe_carn.*(-MEAN_fe_carn+MAX_fe_carn)/(MEAN_fe_carn-MIN_fe_carn));

Spec_NRint_count=zeros(NRspec,1);
delta_list=cell(1,NRspec);
fe_list=cell(1,NRspec);
for IntNR_t=1:NRint_t
    PREY=INT_list_NET_uns(IntNR_t,1);
    PRED=INT_list_NET_uns(IntNR_t,2);
    Spec_NRint_count(PRED,1)=Spec_NRint_count(PRED,1)+1;
    if TL_list(PREY,1)==1 %% herbivore interaction
        delta=rand.*(MAX_delta_herb-MIN_delta_herb)+MIN_delta_herb;
        delta_list{PRED}(1,Spec_NRint_count(PRED,1))=delta;
        
        if rand<FRAC_fe_herb
            fe=MIN_fe_herb+sort(betarnd(a_beta_fe_herb,b_beta_fe_herb,1,1)).*(MAX_fe_herb-MIN_fe_herb);
        else
            fe=1;
        end
        fe_list{PRED}(1,Spec_NRint_count(PRED,1))=fe;
        
    else %% carnivore interaction
        delta=rand.*(MAX_delta_carn-MIN_delta_carn)+MIN_delta_carn;
        delta_list{PRED}(1,Spec_NRint_count(PRED,1))=delta;
        
        if rand<FRAC_fe_carn
            fe=MIN_fe_carn+sort(betarnd(a_beta_fe_carn,b_beta_fe_carn,1,1)).*(MAX_fe_carn-MIN_fe_carn);
        else
            fe=1;
        end
        fe_list{PRED}(1,Spec_NRint_count(PRED,1))=fe;
        
    end
end

%% !!!!!! SORTED STRENGTH LIST !!!!!! - sorting strong interactions with efficient links
INT_ORDER_list=cell(1,NRspec);
DELTA_cvalue=cell(1,NRspec);
DELTA_cvalue_sorted=cell(1,NRspec);
for SpecNR=1:NRspec
    
    DEGREE=Spec_NRprey_count(SpecNR,1);
    if DEGREE>=1
        
        %% random order of assigning relative interaction strengths
        if STRENGTH_SORTING<=2 || STRENGTH_SORTING==4
            INT_ORDER_list{SpecNR}=randperm(DEGREE);
        end
        
        %% strong interactions go together with efficient links
        if STRENGTH_SORTING>=2
            Strth_list{SpecNR}=sort(Strth_list{SpecNR},'descend');
            fe_list{SpecNR}=sort(fe_list{SpecNR}, 'descend');
        end
        
        %% niche order when assigning relative interaction strengths
        if STRENGTH_SORTING==3
            DELTA_cvalue{SpecNR}=abs(NETnichevalue_uns(PREY_list{SpecNR})-NETcvalue_uns(SpecNR,1))';
            DELTA_cvalue_sorted{SpecNR}=sort(DELTA_cvalue{SpecNR});
            
            INT_ORDER_list{SpecNR}=zeros(1,DEGREE);
            for NR_INT_ORDER=1:DEGREE
                POS_NR_INT_ORDER=find(DELTA_cvalue{SpecNR}==DELTA_cvalue_sorted{SpecNR}(1,NR_INT_ORDER));
                INT_ORDER_list{SpecNR}(1,POS_NR_INT_ORDER)=NR_INT_ORDER;
            end
        end
        
        %% strong interactions go together with inefficient links
        if STRENGTH_SORTING==4
            Strth_list{SpecNR}=sort(Strth_list{SpecNR},'descend');
            fe_list{SpecNR}=sort(fe_list{SpecNR}, 'ascend');
        end
        
    end
end

%% Make list of trophic interactions: prey, predator, RELative strength, delta, fe
INT_list_t_REL=INT_list_NET_uns;
Spec_NRint_count=zeros(NRspec,1);
for IntNR_t=1:NRint_t
    
    %% determin INT position on Predator lists
    PRED=INT_list_NET_uns(IntNR_t,2);
    Spec_NRint_count(PRED,1)=Spec_NRint_count(PRED,1)+1;
    INT_POSITION=INT_ORDER_list{PRED}(1,Spec_NRint_count(PRED,1));
    
    INT_list_t_REL(IntNR_t,3)=Strth_list{PRED}(1,INT_POSITION);
    %INT_list_t_REL(IntNR_t,4)=delta_list{PRED}(1,INT_POSITION); %% this should be Spec_NRint_count(PRED,1)?
    INT_list_t_REL(IntNR_t,4)=delta_list{PRED}(1,Spec_NRint_count(PRED,1));
    INT_list_t_REL(IntNR_t,5)=fe_list{PRED}(1,INT_POSITION);
end