function [Neq,Jac,EIGEN,DOM_EIGEN,feasibleCrit,stableCrit]=func_eigenLotVolt(NRspec,RT,A)

%% empty data
Neq=[];
Jac=[];
EIGEN=[];
DOM_EIGEN=[];
stableCrit=0;
feasibleCrit=0;

%% feasTRSH
feasibleTRSH=1e-4;
stableTRSH=0;%-1e-20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% check A is nonsingular %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% clear warnings
lastwarn('')

%% basic Neq calculation
Neq=-(A^-1)*RT;

%% if singular or badly scaled
[~, id_text]=lastwarn;
SINGULAR=strcmp(id_text,'MATLAB:singularMatrix');
NEARSINGULAR=strcmp(id_text,'MATLAB:nearlySingularMatrix');
ILLCOND=strcmp(id_text,'MATLAB:illConditionedMatrix');

MATRIX_COND=0;
if SINGULAR==1 || NEARSINGULAR==1 || ILLCOND==1
    MATRIX_COND=1;
end

if MATRIX_COND==0 && sum(Neq>=feasibleTRSH)==NRspec
    feasibleCrit=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% calculate type I stability %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Jac=diag(RT+sum((A.*(ones(NRspec,1)*Neq')),2))+A.*(Neq*ones(1,NRspec));

if sum(sum(isnan(Jac)))==0 && sum(sum(isinf(Jac)))==0 && MATRIX_COND==0
    EIGEN=eig(Jac);
    SUMEIGEN=sum(EIGEN<=stableTRSH);
    DOM_EIGEN=EIGEN(find(real(EIGEN)==max(real(EIGEN))),1);
    if SUMEIGEN==NRspec
        stableCrit=1;
    end
end