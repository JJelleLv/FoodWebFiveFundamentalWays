function ASSUBSLVFile = analyse_ASSUBS_LV_TPI(NETnr,resultFolder,paramSetName,paramSetFolder,networkFile,printOutput,replaceASSUBSLV)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% FIND Alternative Stable Sub-networks %%%%%%%%%%%%%%%%
%%%%%% Multi-species variety of Yodzis and Innes, 1992 %%%%%
%%%%%% Code by Jelle Lever %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% LV analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load parameter settings
paramSetData=load(sprintf('%s%sPARAMSET_%s',paramSetFolder,filesep,paramSetName));

%% keys
resultsKey_TOPOLOGY=paramSetData.resultsKey_TOPOLOGY;
resultsKey_Neq_M=paramSetData.resultsKey_Neq_M;
resultsKey_INTDIST=paramSetData.resultsKey_INTDIST;

%% load parameters
paramFolder = sprintf('%s%sSTABLE_PARAMETERS_%s%sSTRUCT_PARAMETERS_%s%sINT_DIST_%s%s%d_NICHE_NET_uns', ...
    resultFolder,filesep,resultsKey_TOPOLOGY,filesep,resultsKey_Neq_M,filesep,resultsKey_INTDIST,filesep,NETnr);

networkDATA=load(networkFile);

%%%%%%%%%%%%%%%%%%%%%
%%%% Output File %%%%
%%%%%%%%%%%%%%%%%%%%%

%% determine file and folder to store in
ASSUBSLVFile = sprintf('%s%sASSUBS_LV_TPI', ...
    paramFolder,filesep);

%% skip if this set already exists
if (~replaceASSUBSLV)
    if exist(sprintf('%s.mat',ASSUBSLVFile),'file')==2
        disp('skipped - output loaded from file')
        if printOutput % if NR_ASS>=2
            ASSUBSLVData=load(ASSUBSLVFile);
            fprintf('\nNR alternative Stable States: %d\n',ASSUBSLVData.NR_ASS);
            fprintf('Min. Number of species: %d\n\n',ASSUBSLVData.MIN_SpecNR);
        end
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Assign parameters %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

NRspec=networkDATA.NRspec;
NRBASAL=networkDATA.NRBASAL;

RT=networkDATA.RT;
A=networkDATA.A;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% search for all ASS combinations %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[NR_ASS,EIGEN_found,Neq_found,COMB_found,MIN_SpecNR,MAXAb_Neq_found,SUBSIZE_NRinv_MAT] = code.stability_functions.func_analyse_ASSUBS_LV(resultFolder,NRspec,RT,A);

%% make ASSUBSINFO
ASSUBSINFO=[NETnr NR_ASS MIN_SpecNR MAXAb_Neq_found]; %% to store later

%% display status
if printOutput % if NR_ASS>=2
    fprintf('NR alternative Stable States: %d\n',NR_ASS);
    fprintf('Min. Number of species: %d\n\n',MIN_SpecNR);
end

%% save info found - not omitted
save(ASSUBSLVFile, ...
    'NR_ASS', 'EIGEN_found', 'Neq_found', 'COMB_found', ...
    'MIN_SpecNR', 'MAXAb_Neq_found', 'SUBSIZE_NRinv_MAT','ASSUBSINFO')

