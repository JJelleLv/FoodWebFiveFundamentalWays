function networkFile = makeNETWORK(NETnr,resultFolder,paramSetName,paramSetFolder,replaceNETWORK)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% FIND PARAMETERS FOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Multi-species variety of Yodzis and Innes, 1992 %%%%%
%%%%%% Code by Jelle Lever %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% maximum attempts
MAX_allCrit_ATTNR=1e3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% GENERAL PARAMETER SETTINGS %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load parameters
paramSetData=load(sprintf('%s%sPARAMSET_%s',paramSetFolder,filesep,paramSetName));

%% keys
resultsKey_TOPOLOGY=paramSetData.resultsKey_TOPOLOGY;
resultsKey_Neq_M=paramSetData.resultsKey_Neq_M;
resultsKey_INTDIST=paramSetData.resultsKey_INTDIST;

% Network structure - niche model parameters
S=paramSetData.S;
CONN=paramSetData.CONN; %% to store in networkFile
dCONN=paramSetData.dCONN;
MINNRBASAL=paramSetData.MINNRBASAL; %% to store in networkFile

% target abundance range
MIN_Neq_target=paramSetData.MIN_Neq_target;
MAX_Neq_target=paramSetData.MAX_Neq_target;

%% NRspec
NRspec=paramSetData.NRspec;

%% body-mass ratios
MIN_mRatio_carn=paramSetData.MIN_mRatio_carn;
MAX_mRatio_carn=paramSetData.MAX_mRatio_carn;

%% allometric scaling exponent
b=paramSetData.b;

%% parameters of basal species
MIN_fr=paramSetData.MIN_fr;
MAX_fr=paramSetData.MAX_fr;
MIN_ar=paramSetData.MIN_ar;
MAX_ar=paramSetData.MAX_ar;

%% Competition
MIN_C=paramSetData.MIN_C;
MAX_C=paramSetData.MAX_C;

%% rK ratios
MIN_rKratio=paramSetData.MIN_rKratio;
MAX_rKratio=paramSetData.MAX_rKratio;

%% parameters of consumers/predators (add zeros for basal species)
MIN_fj=paramSetData.MIN_fj;
MAX_fj=paramSetData.MAX_fj;
MIN_aj=paramSetData.MIN_aj;
MAX_aj=paramSetData.MAX_aj;
MIN_at=paramSetData.MIN_at;
MAX_at=paramSetData.MAX_at;

%% distribution of interaction strenghts
FRAC_int_nrnd=paramSetData.FRAC_int_nrnd;
VAR_int=paramSetData.VAR_int;

%% feedinng efficiencies
MIN_delta_herb=paramSetData.MIN_delta_herb;
MAX_delta_herb=paramSetData.MAX_delta_herb;
MIN_delta_carn=paramSetData.MIN_delta_carn;
MAX_delta_carn=paramSetData.MAX_delta_carn;

%% distribution of link efficiency
FRAC_fe_herb=paramSetData.FRAC_fe_herb;
MIN_fe_herb=paramSetData.MIN_fe_herb;
MAX_fe_herb=paramSetData.MAX_fe_herb;
MEAN_fe_herb=paramSetData.MEAN_fe_herb;
VAR_fe_herb=paramSetData.VAR_fe_herb;

FRAC_fe_carn=paramSetData.FRAC_fe_carn;
MIN_fe_carn=paramSetData.MIN_fe_carn;
MAX_fe_carn=paramSetData.MAX_fe_carn;
MEAN_fe_carn=paramSetData.MEAN_fe_carn;
VAR_fe_carn=paramSetData.VAR_fe_carn;

%%%%%%%%%%%%%%%%%%%%
%%%% CONDITIONS %%%%
%%%%%%%%%%%%%%%%%%%%

%% STRENGTH SORTING
STRENGTH_SORTING=paramSetData.STRENGTH_SORTING;

%% minimum average feeding rate
MIN_predSat=paramSetData.MIN_predSat;
MAX_predSat=paramSetData.MAX_predSat;

%% EIGEN value treshold
PARAM_MAX_EIG=paramSetData.PARAM_MAX_EIG;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Output Folder and File %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% determine file and folder to store in
paramFolder = sprintf('%s%sSTABLE_PARAMETERS_%s%sSTRUCT_PARAMETERS_%s%sINT_DIST_%s%s%d_NICHE_NET_uns', ...
    resultFolder,filesep,resultsKey_TOPOLOGY,filesep,resultsKey_Neq_M,filesep,resultsKey_INTDIST,filesep,NETnr);
networkFile = sprintf('%s%sNETWORK_PARAM_TPI', ...
    paramFolder,filesep);

%% skip if this set already exists
if (~replaceNETWORK)
    if exist(sprintf('%s.mat',networkFile),'file')==2        
        disp('skipped - Network already exists')
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Make parameter set %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% search for networks that are feasible and stable
allCrit_ATTNR=0;
while true
    
    %% count times a new network is generated
    allCrit_ATTNR=allCrit_ATTNR+1;
    fprintf('\nAttempt NR: %d\n', allCrit_ATTNR)
    if allCrit_ATTNR>=MAX_allCrit_ATTNR
        error('FAILURE TO GENERATE FEASIBLE/STABLE NETWORK')
    end
    
    %% empty output
    mr=[];
    mc=[];
    Ri=[];
    Ki=[];
    Jk=[];
    Tk=[];
    RT=[];
    A=[];
    Jac=[];
    MATRIX_COND=[];
    EIGEN=[];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Generate NETWORK STRUCTURE %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %disp('> generate network topology')
    
    %% use niche model
    [NETnichevalue_uns,NETnichemax_uns,NETnichemin_uns,NETcvalue_uns,NETrvalue_uns,NET_uns, ...
    NRBASAL,MINNRBASAL,INT_list_NET_uns,BASAL_list,TOPPRED_list,TL_list,NRint_uns,bniche, ...
    CIRCLE_Int_list,CIRCLE_Int_list_NR,CIRCLENR,CIRCLE_nodes,CIRCLE_links,Spec_NRprey_count] = code.generate_network_topology.func_niche_model_MINBS_CIRCLES(S,CONN,dCONN,MINNRBASAL);

    %% Collect basic info about the network
    BASAL_SpecNRs=find(TL_list==1);
    NONBASAL_SpecNRs=find(TL_list>=2);

    NRNONBASAL=length(NONBASAL_SpecNRs);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Sample parameters from distributions %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %disp('> assign parameters to a given network topology')

    %% assign target abundances
    Neq_target=rand(S,1).*(MAX_Neq_target-MIN_Neq_target)+MIN_Neq_target;
    
    %% parameters of basal species
    fr=(TL_list==1).*(rand(NRspec,1).*(MAX_fr-MIN_fr)+MIN_fr);
    ar=(TL_list==1).*(rand(NRspec,1).*(MAX_ar-MIN_ar)+MIN_ar);
    
    %% parameters of non-basal species
    fj=(TL_list>1).*(rand(NRspec,1).*(MAX_fj-MIN_fj)+MIN_fj);
    aj=(TL_list>1).*(rand(NRspec,1).*(MAX_aj-MIN_aj)+MIN_aj);
    at=(TL_list>1).*(rand(NRspec,1).*(MAX_at-MIN_at)+MIN_at);
    
    %% ratio's
    rKratio=(rand(NRspec,1).*(MAX_rKratio-MIN_rKratio)+MIN_rKratio).*(TL_list==1);
    rKratio(rKratio==0)=NaN;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% trophic interactions list - relative feeding preferences %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    NRint_t=NRint_uns;
    
    [INT_list_t_REL]=code.assign_parameters.func_relInt_t(NRspec,NRint_t,NET_uns,INT_list_NET_uns, ...
        NETnichevalue_uns,NETcvalue_uns,TL_list,Spec_NRprey_count, ...
        FRAC_int_nrnd,VAR_int,FRAC_fe_herb,MIN_fe_herb,MAX_fe_herb,MEAN_fe_herb,VAR_fe_herb, ...
        FRAC_fe_carn,MIN_fe_carn,MAX_fe_carn,MEAN_fe_carn,VAR_fe_carn, ...
        MIN_delta_herb,MAX_delta_herb,MIN_delta_carn,MAX_delta_carn,STRENGTH_SORTING);
    
    %% relative feeding matrix - needed for predSat, thetaNeq_mat*totFeedRate-Tk=0
    thetaNeq_mat=zeros(NRspec,NRspec);
    thetaNeq=zeros(NRspec,1);
    for IntNR_t=1:NRint_t
        PREY=INT_list_t_REL(IntNR_t,1);
        PRED=INT_list_t_REL(IntNR_t,2);
        
        thetaNeq(PRED,1)=thetaNeq(PRED,1)+INT_list_t_REL(IntNR_t,3).*Neq_target(PREY,1);
        thetaNeq_mat(PRED,PRED)=thetaNeq_mat(PRED,PRED)+INT_list_t_REL(IntNR_t,3).*Neq_target(PREY,1);
        thetaNeq_mat(PREY,PRED)=thetaNeq_mat(PREY,PRED)-((INT_list_t_REL(IntNR_t,3).*Neq_target(PRED,1))./((1-INT_list_t_REL(IntNR_t,4)).*INT_list_t_REL(IntNR_t,5)));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% competitive interactions list - relative interaction strengths %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [INT_list_c,NRint_c]=code.assign_parameters.func_relInt_c(NRspec,NRBASAL,TL_list,MIN_C,MAX_C);

    %% make kappa matrix (relative competition) - needed for Ri Ki calculation
    kappa=zeros(NRspec,NRspec);
    for IntNR_c=1:NRint_c
        Spec1=INT_list_c(IntNR_c,1);
        Spec2=INT_list_c(IntNR_c,2);
        
        kappa(Spec1,Spec2)=-INT_list_c(IntNR_c,3); %% a minus was added here!
        kappa(Spec2,Spec1)=-INT_list_c(IntNR_c,3);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% assign initial mc body masses %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% initial mc
    mcmin_init=1;
    mcmax_init=100;
    mc=(TL_list>1).*(mcmin_init+rand(NRspec,1)*(mcmax_init-mcmin_init));
    
    %% update mc until predSat within range
    %disp('> update mc until predSat within range')
    [mc,predSatFound]=code.assign_parameters.func_updateMcPredSat(mc,NRspec,NONBASAL_SpecNRs,thetaNeq_mat,thetaNeq,fj,aj,at,b,MIN_predSat,MAX_predSat);
    
    %% update mc until stable
    stableFound=0;
    if predSatFound==1
        %disp('> update mc until stable')
        [mc,stableFound]=code.assign_parameters.func_updateMcStable(mc,NRspec,TL_list,NRBASAL,NONBASAL_SpecNRs,NRint_t,INT_list_t_REL,thetaNeq_mat,thetaNeq,kappa,Neq_target,fj,aj,at,fr,ar,b,rKratio,MIN_predSat,MAX_predSat,PARAM_MAX_EIG);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Determine predator-prey body-mass ratios %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% update mc untill mcRatio ok
    mRatioFound=0;
    if predSatFound==1 && stableFound==1
        %disp('> update mc until mRatio within range')
        [mc,mRatioFound]=code.assign_parameters.func_updateMRatio(mc,NRspec,TL_list,NRBASAL,NONBASAL_SpecNRs,NRint_t,CIRCLE_Int_list,INT_list_t_REL,thetaNeq_mat,thetaNeq,kappa,Neq_target,fj,aj,at,fr,ar,b,rKratio,MIN_predSat,MAX_predSat,PARAM_MAX_EIG,MIN_mRatio_carn,MAX_mRatio_carn);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% END when all criteria fulfilled %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf('predSatFound = %d\n', predSatFound)
    fprintf('stableFound = %d\n', stableFound)
    fprintf('mRatioFound = %d\n', mRatioFound)
    if predSatFound==1 && stableFound==1 && mRatioFound==1
        break
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Check and create output %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% determine predSat, Jk, and Tk
[predSat,predSat_nonbs,totFeedRate,Jk,Tk]=code.assign_parameters.func_determinePredSat(NRspec,mc,fj,aj,at,b,thetaNeq_mat,thetaNeq,NONBASAL_SpecNRs);

%% determine mr, Ri and Ki, and LotVolt
[mr,Ri,Ki,RT,A,Jac,~,EIGEN,stableCrit]=code.assign_parameters.func_assignMrRiKiLotVolt(NRspec,TL_list,NRBASAL,NRint_t,INT_list_t_REL,kappa,totFeedRate,Tk,Neq_target,fr,ar,b,rKratio,PARAM_MAX_EIG);

%% total predSat
dPredSat_nonbs=(predSat_nonbs>MAX_predSat).*(predSat_nonbs-MAX_predSat)+(predSat_nonbs<MIN_predSat).*(MIN_predSat-predSat_nonbs);
sum_dPredSat_nonbs=sum(dPredSat_nonbs);

%% check predSat
if sum_dPredSat_nonbs>0
    error('predSat not in range')
end

%% check stability
if stableCrit==0
    error('not stable')
end

%%%%%%%%%%%%%%%%%%%
%%%% Save DATA %%%%
%%%%%%%%%%%%%%%%%%%

if (~exist(paramFolder, 'dir'))
    mkdir(paramFolder);
end

save(networkFile,'NRspec', ...
    'NETnichevalue_uns','NETnichemax_uns','NETnichemin_uns','NETcvalue_uns','NETrvalue_uns','NET_uns', ...
    'NRBASAL','MINNRBASAL','INT_list_NET_uns','BASAL_list','TOPPRED_list','TL_list','NRint_uns','bniche', ...
    'CIRCLE_Int_list','CIRCLE_Int_list_NR','CIRCLENR','CIRCLE_nodes','CIRCLE_links','Spec_NRprey_count', ...
    'BASAL_SpecNRs','NONBASAL_SpecNRs','NRNONBASAL', ...
    'Neq_target','fr','ar','fj','aj','at','rKratio', ...
    'NRint_t','INT_list_t_REL','thetaNeq_mat','thetaNeq','NRint_c','INT_list_c','kappa', ...
    'mr','Ri','Ki','mc','Jk','Tk','predSat','totFeedRate','RT','A','Jac','EIGEN');

%% display output
fprintf('%d - SUCCESS\n', NETnr)

