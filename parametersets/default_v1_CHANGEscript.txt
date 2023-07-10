function changeSetFile = STEP0_envChangeSettings

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% store parameter settings and generate key to save under %%%%
%%%%%% Code by Jelle Lever, 2016 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SETNAME='default_v1';
replaceFile=false;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CHANGE SETTINGS %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

%% only herb or carnivore interactions
ONLY_herb_carn=0; %% 0=all, 1=only herb, 2=only carn

%% distribution of interaction strenghts
CHANGE_int=1;
NEW_FRAC_int_nrnd=0.1;
NEW_VAR_int=0.0675; % variance for species with degree=2 %% 0.0675

%% minimum average feeding rate
CHANGE_totFeedRate=0;
NEW_RELtotFeedRate_min=0.85;
NEW_RELtotFeedRate_max=1.15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% info about uniform distribution %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% the variance of species with degree=2 when dist is uniform
VAR_int_unif=(1./12).*((1-NEW_FRAC_int_nrnd).^2);

disp(['For a uniform distribution:'])
disp(['VAR_int_unif=', num2str(VAR_int_unif)])
disp([' '])

%%%%%%%%%%%%%%%%%%%
%%%% MAKE KEYS %%%%
%%%%%%%%%%%%%%%%%%%

resultsKey_CHANGE=SETNAME;

% %%% parameter based
% resultsKey_CHANGE=sprintf('%d_%d_%d_%d_%d_%d_%d',ONLY_herb_carn,CHANGE_int,round(NEW_FRAC_int_nrnd*100),round(NEW_VAR_int*1000), ...
%     CHANGE_totFeedRate,round(NEW_RELtotFeedRate_min*100),round(NEW_RELtotFeedRate_max*100))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% STORE PARAMETERS AND KEYS %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

paramSetFolder = sprintf('parametersets');
if (~exist(paramSetFolder, 'dir'))
    mkdir(paramSetFolder);
end

changeSetFile = sprintf('%s%sCHANGESET_%s',paramSetFolder,filesep,SETNAME);

%% save parameters and keys
if exist(sprintf('%s.mat',changeSetFile),'file')==2 && replaceFile==0
    error('CHANGESET already stored under this name/number')
end

save(changeSetFile);
copyfile(['STEP0_envChangeSettings.m'],sprintf('%s%s%s_CHANGEscript.txt', paramSetFolder,filesep,SETNAME));