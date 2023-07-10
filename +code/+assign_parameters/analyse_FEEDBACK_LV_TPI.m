function FEEDBACKLVFile = analyse_FEEDBACK_LV_TPI(NETnr,resultFolder,paramSetName,paramSetFolder,networkFile,printOutput,replaceFEEDBACKLV)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Deteremine total feedback at different levels %%%%%%%
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
FEEDBACKLVFile = sprintf('%s%sFEEDBACK_LV_TPI', ...
    paramFolder,filesep);

%% skip if this set already exists
if (~replaceFEEDBACKLV)
    if exist(sprintf('%s.mat',FEEDBACKLVFile),'file')==2
        disp('skipped - output loaded from file')
        
        % print output
        if printOutput
            FEEDBACKLVData=load(FEEDBACKLVFile);
            fprintf('\nFk (total feedback): [%s]\n', join(string(FEEDBACKLVData.Fk), ','));
            fprintf('Hk (Hurwitz determinants): [%s]\n\n', join(string(FEEDBACKLVData.Hk), ','));
        end
        
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Assign parameters %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

Jac=networkDATA.Jac;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% determine and analyse total feedback %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[Fk]=code.stability_functions.func_feedback(Jac);
[Fk,FkPos,FkNeg]=code.stability_functions.func_feedback_posneg(Jac);
[Hk]=code.stability_functions.func_HurDet(Fk);

% print output
if printOutput
    fprintf('Fk (total feedback): [%s]\n', join(string(Fk), ','));
    fprintf('Hk (Hurwitz determinants): [%s]\n\n', join(string(Hk), ','));
end

%% save info found - not omitted
save(FEEDBACKLVFile,'Fk','FkPos','FkNeg','Hk')

