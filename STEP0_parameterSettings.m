
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Store parameter settings and generate keys to save under %%%%
%%%%%% Code by Jelle Lever %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SETNAME='default_v1';
replaceFile=false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PARAMETER SETTINGS %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Network structure - niche model parameters
S=22;
CONN=0.16; % C=L/S^2, L=number actual links
dCONN=0.02;
MINNRBASAL=1;

%% NRspec
NRspec=S;

% target abundance range
MIN_Neq_target=1.5;
MAX_Neq_target=2.5;

%% body-mass ratios
MIN_mRatio_carn=0.5;
MAX_mRatio_carn=20;%100;

%% allometric scaling exponent
b=-0.25;

%% parameters of basal species
MIN_fr=1;
MAX_fr=1;
MIN_ar=1;
MAX_ar=1;

%% Competition
MIN_C=0.1;
MAX_C=0.5;

%% rK ratios
MIN_rKratio=0.2;
MAX_rKratio=1;

%% parameters of consumers/predators
MIN_fj=1;
MAX_fj=1;
MIN_aj=2.512; %% 2.512 for invertebrates, 3.528 for ectoterm vertebrates
MAX_aj=2.512;
MIN_at=0.314; %% 0.314 for invertebrates, 0.88 for ectoterm vertebrates
MAX_at=0.314;

%% distribution of interaction strenghts
FRAC_int_nrnd=0.1;
VAR_int=0.03; % variance for species with degree=2 %% 0.0675

%% feedinng efficiencies
MIN_delta_herb=0.4;%0; %0.55; %% 0.55 for herbivore interactions (high means inefficient)
MAX_delta_herb=0.7;%MIN_delta_herb; %0.55;
MIN_delta_carn=0.05;%MIN_delta_herb; %0.15; %% 0.15 for carnivore interactions
MAX_delta_carn=0.25;%MIN_delta_herb; %0.15;

%% distribution of link efficiency
FRAC_fe_herb=1; %% the raction of links that are inefficient (fe~=1)
MIN_fe_herb=0.1; %0.2
MAX_fe_herb=1; %1, max is 1
MEAN_fe_herb=0.75; % 0.6
VAR_fe_herb=0.11; %0.053333 %<<<<<

FRAC_fe_carn=FRAC_fe_herb;
MIN_fe_carn=MIN_fe_herb; % 0.2
MAX_fe_carn=MAX_fe_herb; % 1, max is 1
MEAN_fe_carn=MEAN_fe_herb; % 0.6
VAR_fe_carn=VAR_fe_herb; %0.053333

%% STRENGTH SORTING
STRENGTH_SORTING=3; %% 1=random, 2=strong int are efficient, 3=niche order

%% minimum average feeding rate
MIN_predSat=0.05;
MAX_predSat=0.75;

%% EIGEN value treshold, slightly smaller than zero (!)
PARAM_MAX_EIG=-1e-4;%-1e-3; %% not smaller than -1e-7

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% info about uniform distribution %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% the variance of species with degree=2 when dist is uniform
VAR_int_unif=(1./12).*((1-FRAC_int_nrnd).^2);

%% the variance and that belongs to a uniform distribution
MEAN_fe_unif_herb=MIN_fe_herb+0.5.*(MAX_fe_herb-MIN_fe_herb);
VAR_fe_unif_herb=1/(((2).^2).*(2+1)).*((MAX_fe_herb-MIN_fe_herb)^2);

MEAN_fe_unif_carn=MIN_fe_carn+0.5.*(MAX_fe_carn-MIN_fe_carn);
VAR_fe_unif_carn=1/(((2).^2).*(2+1)).*((MAX_fe_carn-MIN_fe_carn)^2);

disp(['For a uniform distribution:'])
disp(['VAR_int_unif=', num2str(VAR_int_unif)])
disp(['MEAN_fe_unif_herb=', num2str(MEAN_fe_unif_herb)])
disp(['VAR_fe_unif_herb=', num2str(VAR_fe_unif_herb)])
disp(['MEAN_fe_unif_carn=', num2str(MEAN_fe_unif_carn)])
disp(['VAR_fe_unif_carn=', num2str(VAR_fe_unif_carn)])
disp([' '])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% make a database of all species combinations %%%%
%%%% to considerably speed up further analysis %%%%%%
%%%% (may produce very large files) %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% result folder
resultFolder = 'results';

fprintf('\nSTORE ALL COMBINATIONS\n')
code.all_combns.all_combns(NRspec,resultFolder);

%%%%%%%%%%%%%%%%%%%
%%%% MAKE KEYS %%%%
%%%%%%%%%%%%%%%%%%%

%%% name based
resultsKey_TOPOLOGY=SETNAME;

resultsKey_Neq_M=SETNAME;

resultsKey_INTDIST=SETNAME;

%%% parameter based
% resultsKey_TOPOLOGY=sprintf('%d_%d_%d_%d',S,round(CONN*100),round(dCONN*100),MINNRBASAL)
% 
% resultsKey_Neq_M=sprintf('%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d',NRspec,round(MIN_Neq_target*100),round(MAX_Neq_target*100), ...
%     round(MIN_mRatio_carn*1e3),round(MAX_mRatio_carn*1e3), ...
%     round(MIN_fr*100),round(MAX_fr*100),round(MIN_ar*100),round(MAX_ar*100), ...
%     round(MIN_fj*100),round(MAX_fj*100),round(MIN_aj*100),round(MAX_aj*100), ...
%     round(MIN_at*100),round(MAX_at*100),round(-PARAM_MAX_EIG*1e7))
% 
% resultsKey_INTDIST=sprintf('%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d',STRENGTH_SORTING,round(FRAC_int_nrnd*100),round(VAR_int*1000), ...
%     round(MIN_delta_herb*100),round(MAX_delta_herb*100),round(MIN_delta_carn*100),round(MAX_delta_carn*100), ...
%     round(FRAC_fe_herb*100),round(MIN_fe_herb*100),round(MAX_fe_herb*100),round(MEAN_fe_herb*100),round(VAR_fe_herb*1000), ...
%     round(FRAC_fe_carn*100),round(MIN_fe_carn*100),round(MAX_fe_carn*100),round(MEAN_fe_carn*100),round(VAR_fe_carn*1000), ...
%     round(MIN_predSat*100),round(MAX_predSat*100), ...
%     round(MIN_C*100),round(MAX_C*100),round(MIN_rKratio*1000),round(MAX_rKratio*1000))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% STORE PARAMETERS AND KEYS %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

paramSetFolder = sprintf('parametersets');
if (~exist(paramSetFolder, 'dir'))
    mkdir(paramSetFolder);
end

paramSetFile = sprintf('%s%sPARAMSET_%s',paramSetFolder,filesep,SETNAME);

%% save parameters and keys
if exist(sprintf('%s.mat',paramSetFile),'file')==2 && replaceFile==0
    error('PARAMSET already stored under this name/number')
end

save(paramSetFile);
copyfile(['STEP0_parameterSettings.m'],sprintf('%s%s%s_PARAMscript.txt', paramSetFolder,filesep,SETNAME));