function [mr,Ri,Ki,RT,A,Jac,MATRIX_COND,EIGEN,stableCrit]=func_assignMrRiKiLotVolt(NRspec,TL_list,NRBASAL,NRint_t,INT_list_t_REL,kappa,totFeedRate,Tk,Neq_target,fr,ar,b,rKratio,PARAM_MAX_EIG)

%% empty data
EIGEN=[];
stableCrit=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% determine mr, Ri and Ki %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off', 'MATLAB:singularMatrix')
warning('off', 'MATLAB:nearlySingularMatrix')
warning('off', 'MATLAB:illConditionedMatrix')

%% make feeding matrix
A_t=zeros(NRspec,NRspec);
for IntNR_t=1:NRint_t
    PREY=INT_list_t_REL(IntNR_t,1);
    PRED=INT_list_t_REL(IntNR_t,2);
    A_t(PRED,PREY)=A_t(PRED,PREY)+totFeedRate(PRED,1).*INT_list_t_REL(IntNR_t,3);
    A_t(PREY,PRED)=A_t(PREY,PRED)-(totFeedRate(PRED,1).*INT_list_t_REL(IntNR_t,3))./((1-INT_list_t_REL(IntNR_t,4)).*INT_list_t_REL(IntNR_t,5));
end

%% calculate feeding on basal
FEEDING_on_BASAL=(TL_list==1).*(A_t*Neq_target);
FEEDING_on_BASAL(TL_list>=2)=NaN;

%% if COMB_FEEDING_on_BASAL are positive, K cannot be correct (=negative).
if sum(FEEDING_on_BASAL<=0)~=NRBASAL
    disp('something wrong with FEEDING on BASAL')
    FEEDING_on_BASAL
    mr=[];
    Ri=[];
    Ki=[];
    RT=[];
    A=[];
    Jac=[];
    MATRIX_COND=1;
    EIGEN=[];
    stableCrit=0;
    return
    error('something wrong with FEEDING on BASAL')
end

%% relation between FEEDING_on_BASAL and R and K
Ri=(TL_list==1).*(-(kappa*Neq_target).*rKratio-FEEDING_on_BASAL);
Ki=(TL_list==1).*(Ri./rKratio);
Ri(isnan(Ri))=0;
Ki(isnan(Ki))=0;

if (sum(Ri<0)+sum(Ki<0))>=1    
    disp('something wrong, negative R and K')
    FEEDING_on_BASAL
    mr=[];
    Ri=[];
    Ki=[];
    RT=[];
    A=[];
    Jac=[];
    MATRIX_COND=1;
    EIGEN=[];
    stableCrit=0;
    return
    error('negative R and K')
end

%% determine prey body mass - Ri=fr.*ar.*(mr.^b)
mr=(Ri./(fr.*ar)).^(1./b);
mr(isnan(mr))=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% make final RT and A %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% growth rates
RT=zeros(NRspec,1);
RT(isnan(Ri)==0)=Ri(isnan(Ri)==0);
RT(isnan(Tk)==0)=-Tk(isnan(Tk)==0);

%% interaction strengths
A_c=kappa.*((Ri./Ki)*ones(1,NRspec));
A_c(isnan(A_c))=0;
A=A_c+A_t;

%% check that there is no change at the eq
DIFF=RT+A*Neq_target;
if sum(abs(DIFF)>=1e-7)>=1
    disp('something is wrong!!!! DIFF~=0')
    mr=[];
    Ri=[];
    Ki=[];
    RT=[];
    A=[];
    Jac=[];
    MATRIX_COND=1;
    EIGEN=[];
    stableCrit=0;
    return
    error('something is wrong!!!! DIFF~=0')
end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% calculate type I stability %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Jac=diag(RT+sum((A.*(ones(NRspec,1)*Neq')),2))+A.*(Neq*ones(1,NRspec));

if sum(sum(isnan(Jac)))==0 && sum(sum(isinf(Jac)))==0 && MATRIX_COND==0
    EIGEN=eig(Jac);
    SUMEIGEN=sum(EIGEN<=PARAM_MAX_EIG);
    if SUMEIGEN==NRspec
        stableCrit=1;
    end
end