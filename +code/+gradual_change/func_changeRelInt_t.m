function [INT_list_t_REL]=func_changeRelInt_t(NRspec,NRint_t,INT_list_t_REL_null,TL_list,NEW_FRAC_int_nrnd,NEW_VAR_int,ONLY_herb_carn)

INT_list_t_REL=INT_list_t_REL_null;

%% determine (sum of) interctions that should be changed    
changeIntList=zeros(NRint_t,1);
changeIntSum=zeros(NRspec,1);
changeSpec_NRprey_count=zeros(NRspec,1);
if ONLY_herb_carn==0
    changeIntList=ones(NRint_t,1);
    changeIntSum=ones(NRspec,1).*(TL_list>=2); 
    for IntNR_t=1:NRint_t
        PRED=INT_list_t_REL_null(IntNR_t,2);
        changeSpec_NRprey_count(PRED,1)=changeSpec_NRprey_count(PRED,1)+1;
    end
elseif ONLY_herb_carn==1
    for IntNR_t=1:NRint_t
        PREY=INT_list_t_REL_null(IntNR_t,1);
        PRED=INT_list_t_REL_null(IntNR_t,2);
        if TL_list(PREY,1)==1
            changeIntList(IntNR_t,1)=1;
            changeIntSum(PRED,1)=changeIntSum(PRED,1)+INT_list_t_REL_null(IntNR_t,3);
            changeSpec_NRprey_count(PRED,1)=changeSpec_NRprey_count(PRED,1)+1;
        end
    end
elseif ONLY_herb_carn==2
    for IntNR_t=1:NRint_t
        PREY=INT_list_t_REL_null(IntNR_t,1);
        PRED=INT_list_t_REL_null(IntNR_t,2);
        if TL_list(PREY,1)>=2
            changeIntList(IntNR_t,1)=1;
            changeIntSum(PRED,1)=changeIntSum(PRED,1)+INT_list_t_REL_null(IntNR_t,3);
            changeSpec_NRprey_count(PRED,1)=changeSpec_NRprey_count(PRED,1)+1;
        end
    end
end

% changeIntList
% changeIntSum
% error('till lalalalala')

%% take relative interaction strengths from dirichlet distribution
A_gamma=-(-1+4.*NEW_VAR_int+2.*NEW_FRAC_int_nrnd-(NEW_FRAC_int_nrnd.^2))./(8.*NEW_VAR_int);

Strth_list=cell(1,NRspec);
for SpecNR=1:NRspec
    
    DEGREE=changeSpec_NRprey_count(SpecNR,1);
    if DEGREE>=1
        
        MIN_int=NEW_FRAC_int_nrnd./DEGREE;
        a_gamma_MAT=ones(1,DEGREE).*A_gamma;
        
        DIST_DIRI=gamrnd(a_gamma_MAT,1,1,DEGREE);
        DIST_DIRI=DIST_DIRI./sum(DIST_DIRI);
        DIST_DIRI=MIN_int+DIST_DIRI.*((1-NEW_FRAC_int_nrnd));
        
        Strth_list{SpecNR}=DIST_DIRI;
        
    end
    
end

%% !!!!!! SORTED STRENGTH LIST !!!!!! - sorting strong interactions with efficient links
INT_ORDER_list=cell(1,NRspec);
for SpecNR=1:NRspec
    
    DEGREE=changeSpec_NRprey_count(SpecNR,1);
    if DEGREE>=1
        
        %% random order of assigning relative interaction strengths
        INT_ORDER_list{SpecNR}=randperm(DEGREE);
        
        %% other sortings are removed from this function 
        % because makes no sense when only herb or carn?
        
    end
end

%% Make list of trophic interactions: prey, predator, RELative strength, delta, fe
changeSpec_NRint_count=zeros(NRspec,1);
for IntNR_t=1:NRint_t
    
    if changeIntList(IntNR_t,1)==1
        
        %% determin INT position on Predator lists
        PRED=INT_list_t_REL_null(IntNR_t,2);
        changeSpec_NRint_count(PRED,1)=changeSpec_NRint_count(PRED,1)+1;
        INT_POSITION=INT_ORDER_list{PRED}(1,changeSpec_NRint_count(PRED,1));
        
        INT_list_t_REL(IntNR_t,3)=Strth_list{PRED}(1,INT_POSITION).*changeIntSum(PRED,1);
        
    end
end