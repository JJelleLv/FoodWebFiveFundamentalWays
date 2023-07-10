function changeFile = changeNETWORK(NETnr,NRchange,NRsteps_full,resultFolder,paramSetName,paramSetFolder,changeSetName,networkFile,replaceChangeNETWORK)

%% load parameter settings
paramSetData=load(sprintf('%s%sPARAMSET_%s',paramSetFolder,filesep,paramSetName));
changeSetData=load(sprintf('%s%sCHANGESET_%s', paramSetFolder,filesep,changeSetName));

%% keys
resultsKey_TOPOLOGY=paramSetData.resultsKey_TOPOLOGY;
resultsKey_Neq_M=paramSetData.resultsKey_Neq_M;
resultsKey_INTDIST=paramSetData.resultsKey_INTDIST;
resultsKey_CHANGE=changeSetData.resultsKey_CHANGE;

%% load parameters
paramFolder = sprintf('%s%sSTABLE_PARAMETERS_%s%sSTRUCT_PARAMETERS_%s%sINT_DIST_%s%s%d_NICHE_NET_uns', ...
    resultFolder,filesep,resultsKey_TOPOLOGY,filesep,resultsKey_Neq_M,filesep,resultsKey_INTDIST,filesep,NETnr);
changeFolder = sprintf('%s%sCHANGESET_%s',paramFolder,filesep,resultsKey_CHANGE);

networkDATA=load(networkFile);

%%%%%%%%%%%%%%%%%%%%%
%%%% Output File %%%%
%%%%%%%%%%%%%%%%%%%%%

%% determine file and folder to store in
changeFile = sprintf('%s%sCHANGE_LV_TPI_%d_%d', ...
    changeFolder,filesep,NRchange,NRsteps_full);

%% skip if this set already exists
if (~replaceChangeNETWORK)
    if exist(sprintf('%s.mat',changeFile),'file')==2
        disp('skipped -  output loaded from file')
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Assign parameters %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

NRspec=networkDATA.NRspec;
TL_list=networkDATA.TL_list;
Spec_NRprey_count=networkDATA.Spec_NRprey_count;

mr=networkDATA.mr;
mc=networkDATA.mc;
Ri=networkDATA.Ri;
Ki=networkDATA.Ki;
Tk=networkDATA.Tk;
NRint_t=networkDATA.NRint_t;
NRint_c=networkDATA.NRint_c;
INT_list_c=networkDATA.INT_list_c;

%% null network
INT_list_t_REL_null=networkDATA.INT_list_t_REL;
totFeedRate_null=networkDATA.totFeedRate;
RT_null=networkDATA.RT;
A_null=networkDATA.A;

%% change info
ONLY_herb_carn=changeSetData.ONLY_herb_carn;
CHANGE_int=changeSetData.CHANGE_int;
NEW_FRAC_int_nrnd=changeSetData.NEW_FRAC_int_nrnd;
NEW_VAR_int=changeSetData.NEW_VAR_int;
CHANGE_totFeedRate=changeSetData.CHANGE_totFeedRate;
NEW_RELtotFeedRate_min=changeSetData.NEW_RELtotFeedRate_min;
NEW_RELtotFeedRate_max=changeSetData.NEW_RELtotFeedRate_max;

%% empty data
INT_list_t_REL_end=cell(1,NRchange);
totFeedRate_end=cell(1,NRchange);
RT_end=cell(1,NRchange);
A_end=cell(1,NRchange);
RT_change=cell(1,NRchange);
A_change=cell(1,NRchange);
Tpp_found=NaN(NRchange,1);
Tpp_stepNR=NaN(NRchange,1);
Tpp_Neq=cell(1,NRchange);
Tpp_minNeq=NaN(NRchange,1);
Tpp_Jac=cell(1,NRchange);
Tpp_EIGEN=cell(1,NRchange);
Tpp_DOM_EIGEN=cell(1,NRchange);
Tpp_feasibleCrit=NaN(NRchange,1);
Tpp_stableCrit=NaN(NRchange,1);
Tpp_HOPF=NaN(NRchange,1);
Tpp_L0=NaN(NRchange,1);
Tpp_L0_MATCONT=NaN(NRchange,1);
NTReq_series=cell(1,NRchange);
NTR_EIGEN_series=cell(1,NRchange);
NTR_DOM_EIGEN_series=cell(1,NRchange);
NRsteps=NaN(NRchange,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% apply global change %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for changeNR=1:NRchange
    
    %% apply change
    [INT_list_t_REL_end{changeNR},totFeedRate_end{changeNR},RT_end{changeNR},A_end{changeNR},RT_change{changeNR},A_change{changeNR},Tpp_found(changeNR,1),Tpp_stepNR(changeNR,1),Tpp_Neq{changeNR},Tpp_minNeq(changeNR,1),Tpp_Jac{changeNR},Tpp_EIGEN{changeNR},Tpp_DOM_EIGEN{changeNR}, ...
        Tpp_feasibleCrit(changeNR,1),Tpp_stableCrit(changeNR,1),Tpp_HOPF(changeNR,1),Tpp_L0(changeNR,1),Tpp_L0_MATCONT(changeNR,1),NTReq_series{changeNR},NTR_EIGEN_series{changeNR},NTR_DOM_EIGEN_series{changeNR},NRsteps(changeNR,1)]=code.gradual_change.func_applyChange(NETnr,changeNR, ...
        NRsteps_full,NRspec,TL_list,Spec_NRprey_count,mr,mc,Ri,Ki,Tk,NRint_t,NRint_c,INT_list_c,INT_list_t_REL_null,totFeedRate_null,RT_null,A_null,ONLY_herb_carn,CHANGE_int,NEW_FRAC_int_nrnd,NEW_VAR_int,CHANGE_totFeedRate,NEW_RELtotFeedRate_min,NEW_RELtotFeedRate_max);

end

%%%%%%%%%%%%%%%%%%%
%%%% Save DATA %%%%
%%%%%%%%%%%%%%%%%%%

if (~exist(changeFolder, 'dir'))
    mkdir(changeFolder);
end

%save(changeFile,'ONLY_herb_carn','CHANGE_int','NEW_FRAC_int_nrnd','NEW_VAR_int','CHANGE_totFeedRate','NEW_RELtotFeedRate_min','NEW_RELtotFeedRate_max','NRsteps_full', ...
%    'INT_list_t_REL_end','totFeedRate_end','RT_end','A_end','RT_change','A_change','Tpp_found','Tpp_stepNR','Tpp_Neq','Tpp_Jac','Tpp_EIGEN','Tpp_DOM_EIGEN', ...
%    'Tpp_feasibleCrit','Tpp_stableCrit','Tpp_HOPF','Tpp_L0','Tpp_L0_MATCONT','NTReq_series','NTR_EIGEN_series','NTR_DOM_EIGEN_series','NRsteps');

save(changeFile,'ONLY_herb_carn','CHANGE_int','NEW_FRAC_int_nrnd','NEW_VAR_int','CHANGE_totFeedRate','NEW_RELtotFeedRate_min','NEW_RELtotFeedRate_max','NRsteps_full', ...
    'INT_list_t_REL_end','totFeedRate_end','RT_end','A_end','RT_change','A_change','Tpp_found','Tpp_stepNR','Tpp_Neq','Tpp_minNeq','Tpp_Jac','Tpp_EIGEN','Tpp_DOM_EIGEN', ...
    'Tpp_feasibleCrit','Tpp_stableCrit','Tpp_HOPF','Tpp_L0','Tpp_L0_MATCONT'); %% without series to save storage space
