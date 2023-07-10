function MAIN_makeAnalyseFoodWebs(NETnr_min,NETnr_max)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MAIN code to make and analyse foodwebs %%%%%%%%%%%%%
%%%% Niche model (Williams and Martinez, 2000) %%%%%%%%%%
%%%% Multi-species variety of Yodzis and Innes, 1992 %%%%
%%%% Analysis of total feedback, alternative %%%%%%%%%%%%
%%%% stable states, and bifurcation points %%%%%%%%%%%%%%
%%%% Code by Jelle Lever %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng('shuffle');
clc

%%% print output
printOutput=true;

%%% result folder
resultFolder = 'results';

%%% parameter settings
paramSetNameSeries{1}='default_v1';
%% paramSetNameSeries{2}= ... further parameter settings can be added

%%% global environmental change
changeSetNameSeries{1}='default_v1';
%% changeSetNameSeries{2}=... further change scenarios can be added

for paramSetNameNR=1:length(paramSetNameSeries)

    for NETnr=NETnr_min:NETnr_max
        
        % record time
        tic
        
        % load parameter settings
        paramSetName=paramSetNameSeries{paramSetNameNR};
        loadPARAMSET=load(sprintf('parametersets%sPARAMSET_%s',filesep,paramSetName));
        paramSetFolder=loadPARAMSET.paramSetFolder;
        fprintf(['\nParameter settings: ', paramSetName,'\n'])
        
        %%% STEP1: Generate network
        fprintf('\nGENERATE NETWORK (paramSetNR, NETnr): %d, %d\n',paramSetNameNR,NETnr)
        
        replaceNETWORK = false;
        networkFile = code.assign_parameters.makeNETWORK(NETnr,resultFolder,paramSetName,paramSetFolder,replaceNETWORK);
        
        %%% STEP2: ANALYSE INITIAL CONDITION3S - LV
        fprintf('\nANALYSE INITIAL CONDITIONS (paramSetNR, NETnr): %d, %d\n',paramSetNameNR,NETnr)
        
        % determine Feedbacks
        disp('determine FEEDBACKS')
        replaceFEEDBACKLV = false;
        FEEDBACKLVFile = code.assign_parameters.analyse_FEEDBACK_LV_TPI(NETnr,resultFolder,paramSetName,paramSetFolder,networkFile,printOutput,replaceFEEDBACKLV);
        
        % Check alternative stable states
        disp('check ASSUBS LV')
        replaceASSUBSLV = false;
        ASSUBSLVFile = code.assign_parameters.analyse_ASSUBS_LV_TPI(NETnr,resultFolder,paramSetName,paramSetFolder,networkFile,printOutput,replaceASSUBSLV);
        
        %%% STEP3: APPLY GLOBAL CHANGE - LV
        fprintf('\nAPPLY GLOBAL CHANGE (paramSetNR, NETnr): %d, %d\n',paramSetNameNR,NETnr)
        for changeSetNameNR=1:length(changeSetNameSeries)
            
            changeSetName=changeSetNameSeries{changeSetNameNR};
            fprintf(['\nChange settings: ', changeSetName,'\n'])
            
            %%% STEP3a: apply global environmental change, check TPP TYPES
            replaceChangeNETWORK = false;
            NRchange=100;
            NRsteps=10000;
            changeFile = code.gradual_change.changeNETWORK(NETnr,NRchange,NRsteps,resultFolder,paramSetName,paramSetFolder,changeSetName,networkFile,replaceChangeNETWORK);

            %%% STEP 3b: make prescision simulation around Tpp
            replaceChangeNETWORKPrecision=false;
            NRprecisionSteps=100;
            NRchangeASSAnalysis=5;
            changePrecisionFile = code.gradual_change.changeNETWORKPrecision(NETnr,NRchange,NRchangeASSAnalysis,NRsteps,NRprecisionSteps,resultFolder,paramSetName,paramSetFolder,changeSetName,networkFile,printOutput,replaceChangeNETWORKPrecision);
            
        end
        
        % end record time
        toc
        
    end
end