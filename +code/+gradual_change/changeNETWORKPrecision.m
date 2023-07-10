function changePrecisionFile = changeNETWORKPrecision(NETnr,NRchange,NRchangeASSAnalysis,NRsteps_full,NRsteps_fullPrecision,resultFolder,paramSetName,paramSetFolder,changeSetName,networkFile,printOutput,replaceChangeNETWORK)

%% load parameter settings
paramSetData=load(sprintf('%s%sPARAMSET_%s',paramSetFolder,filesep,paramSetName));
changeSetData=load(sprintf('%s%sCHANGESET_%s',paramSetFolder,filesep,changeSetName));

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
changePrecisionFile = sprintf('%s%sCHANGEPrecision_LV_TPI_%d_%d_%d', ...
    changeFolder,filesep,NRchange,NRsteps_full,NRsteps_fullPrecision);

%% skip if this set already exists
if (~replaceChangeNETWORK)
    if exist(sprintf('%s.mat',changePrecisionFile),'file')==2
        disp('skipped - output loaded from file')
        
        if printOutput
            changePrecisionData=load(changePrecisionFile);
            Tpp_found_Precision=changePrecisionData.Tpp_found_Precision;
            Tpp_feasibleCrit=changePrecisionData.Tpp_feasibleCrit;
            Tpp_TRSCRT=changePrecisionData.Tpp_TRSCRT;
            Tpp_HOPF=changePrecisionData.Tpp_HOPF;
            Tpp_L0=changePrecisionData.Tpp_L0;
            Tpp_M=changePrecisionData.Tpp_M;
            Tpp_HOPF_NR_ASS=changePrecisionData.Tpp_HOPF_NR_ASS;
            Tpp_HOPF_MIN_SpecNR=changePrecisionData.Tpp_HOPF_MIN_SpecNR;
            
            for changeNR=1:NRchange
                if Tpp_found_Precision(changeNR,1)==1
                    if Tpp_feasibleCrit(changeNR,1)==0
                        if Tpp_TRSCRT(changeNR,1)==1
                            fprintf('%d - %d - Transcritical bicurcation at E=%.4f\n',NETnr,changeNR,Tpp_M(changeNR,1))
                        else
                            fprintf('%d - %d - Loss of feasibility at E=%.4f\n',NETnr,changeNR,Tpp_M(changeNR,1))
                            %Tpp_Neq
                            %Tpp_EIGEN
                            %stableCrit_TRSCRTSUB
                            %NRinv
                            %pause
                        end
                    else
                        if Tpp_HOPF(changeNR,1)==0
                            fprintf('>>> %d - %d - STRANGE OUTPUT - FEASIBLE - NO HOPF!!!! E=%.4f\n',NETnr,changeNR,Tpp_M(changeNR,1))
                        else
                            if Tpp_L0(changeNR,1)<0
                                fprintf('%d - %d - Supercritical HOPF (limit cycle), L0: %.4f at E=%.4f, NRASS=%d, minSpecNR=%d\n',NETnr,changeNR,Tpp_L0(changeNR,1),Tpp_M(changeNR,1),Tpp_HOPF_NR_ASS(changeNR,1),Tpp_HOPF_MIN_SpecNR(changeNR,1))
                            elseif Tpp_L0(changeNR,1)>0
                                fprintf('%d - %d - Subcritical HOPF (instability), L0: %.4f at E=%.4f, NRASS=%d, minSpecNR=%d\n',NETnr,changeNR,Tpp_L0(changeNR,1),Tpp_M(changeNR,1),Tpp_HOPF_NR_ASS(changeNR,1),Tpp_HOPF_MIN_SpecNR(changeNR,1))
                            else
                                fprintf('%d - %d - Unknown HOPF (likely sum(EIGEN>0)>2), L0: %.4f at E=%.4f, NRASS=%d, minSpecNR=%d\n',NETnr,changeNR,Tpp_L0(changeNR,1),Tpp_M(changeNR,1),Tpp_HOPF_NR_ASS(changeNR,1),Tpp_HOPF_MIN_SpecNR(changeNR,1))
                            end
                        end
                    end
                else
                    fprintf('%d - %d - No Tipping point FOUND!!\n',NETnr,changeNR)
                end
            end
        end
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Load changeFile %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

changeFile = sprintf('%s%sCHANGE_LV_TPI_%d_%d', ...
    changeFolder,filesep,NRchange,NRsteps_full);

changeDATA=load(changeFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Assign parameters %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

NRspec=networkDATA.NRspec;

%% null network
RT_null=networkDATA.RT;
A_null=networkDATA.A;

%% change info
RT_change=changeDATA.RT_change;
A_change=changeDATA.A_change;

%% Tpp info
Tpp_stepNR=changeDATA.Tpp_stepNR;
Tpp_found=changeDATA.Tpp_found;

%% empty data
Tpp_found_Precision=NaN(NRchange,1);
Tpp_stepNR_Precision=NaN(NRchange,1);
Tpp_Neq=cell(1,NRchange);
Tpp_minNeq=NaN(NRchange,1);
Tpp_Jac=cell(1,NRchange);
Tpp_EIGEN=cell(1,NRchange);
Tpp_DOM_EIGEN=cell(1,NRchange);
Tpp_feasibleCrit=NaN(NRchange,1);
Tpp_stableCrit=NaN(NRchange,1);
Tpp_TRSCRT=NaN(NRchange,1);
Tpp_HOPF=NaN(NRchange,1);
Tpp_L0=NaN(NRchange,1);
Tpp_L0_MATCONT=NaN(NRchange,1);
Tpp_HOPF_NR_ASS=NaN(NRchange,1);
Tpp_HOPF_MIN_SpecNR=NaN(NRchange,1);
Tpp_M=NaN(NRchange,1);
NTReq_series=cell(1,NRchange);
NTR_EIGEN_series=cell(1,NRchange);
NTR_DOM_EIGEN_series=cell(1,NRchange);
NRsteps_Precision=NaN(NRchange,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% apply global change %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for changeNR=1:NRchange
    
    checkStepNR=Tpp_stepNR(changeNR,1)-1;
    if Tpp_found(changeNR,1)==0 %% it could be that this is either instability within the first step or no instability at all
        checkStepNR=NRsteps_full;
    end
    
    %% apply change
    [Tpp_found_Precision(changeNR,1),Tpp_stepNR_Precision(changeNR,1),Tpp_Neq{changeNR},Tpp_minNeq(changeNR,1),Tpp_Jac{changeNR},Tpp_EIGEN{changeNR},Tpp_DOM_EIGEN{changeNR}, ...
        Tpp_feasibleCrit(changeNR,1),Tpp_stableCrit(changeNR,1),Tpp_TRSCRT(changeNR,1),Tpp_HOPF(changeNR,1),Tpp_L0(changeNR,1),Tpp_L0_MATCONT(changeNR,1),Tpp_HOPF_NR_ASS(changeNR,1),Tpp_HOPF_MIN_SpecNR(changeNR,1),NTReq_series{changeNR},NTR_EIGEN_series{changeNR},NTR_DOM_EIGEN_series{changeNR},NRsteps_Precision(changeNR,1),Tpp_M(changeNR,1)]=code.gradual_change.func_applyChangePrecision(NETnr,changeNR,NRchangeASSAnalysis,NRsteps_full,NRsteps_fullPrecision,NRspec,RT_null,A_null,RT_change{changeNR},A_change{changeNR},resultFolder,checkStepNR);
        
end

%%%%%%%%%%%%%%%%%%%
%%%% Save DATA %%%%
%%%%%%%%%%%%%%%%%%%

if (~exist(changeFolder, 'dir'))
    mkdir(changeFolder);
end

% save(changePrecisionFile,'NRsteps_full','RT_null','A_null','RT_change','A_change', ...
%     'Tpp_found_Precision','Tpp_stepNR_Precision','Tpp_Neq','Tpp_minNeq','Tpp_Jac','Tpp_EIGEN','Tpp_DOM_EIGEN', ...
%     'Tpp_feasibleCrit','Tpp_stableCrit','Tpp_TRSCRT','Tpp_HOPF','Tpp_L0','Tpp_L0_MATCONT', ...
%     'NTReq_series','NTR_EIGEN_series','NTR_DOM_EIGEN_series','NRsteps_Precision','Tpp_M');

save(changePrecisionFile,'NRchangeASSAnalysis','NRsteps_full','RT_null','A_null','RT_change','A_change', ...
    'Tpp_found_Precision','Tpp_stepNR_Precision','Tpp_Neq','Tpp_minNeq','Tpp_Jac','Tpp_EIGEN','Tpp_DOM_EIGEN', ...
    'Tpp_feasibleCrit','Tpp_stableCrit','Tpp_TRSCRT','Tpp_HOPF','Tpp_L0','Tpp_L0_MATCONT','Tpp_HOPF_NR_ASS','Tpp_HOPF_MIN_SpecNR', ...
    'NRsteps_Precision','Tpp_M'); %% without series to save storage space