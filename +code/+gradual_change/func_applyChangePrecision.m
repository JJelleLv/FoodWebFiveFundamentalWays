function [Tpp_found,Tpp_stepNR_Precision,Tpp_Neq,Tpp_minNeq,Tpp_Jac,Tpp_EIGEN,Tpp_DOM_EIGEN, ...
Tpp_feasibleCrit,Tpp_stableCrit,Tpp_TRSCRT,Tpp_HOPF,Tpp_L0,Tpp_L0_MATCONT,Tpp_HOPF_NR_ASS,Tpp_HOPF_MIN_SpecNR,NTReq_series,NTR_EIGEN_series,NTR_DOM_EIGEN_series,NRsteps_fullPrecision,Tpp_M]=func_applyChangePrecision(NETnr,changeNR,NRchangeASSAnalysis,NRsteps_full,NRsteps_fullPrecision,NRspec,RT_null,A_null,RT_change,A_change,resultFolder,stepNR)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MAKE EMPTY DATA SERIES %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% stability properties
NTReq_series=[];
NTR_EIGEN_series=[];
NTR_DOM_EIGEN_series=[];

%% first unstable point
Tpp_Neq=[];
Tpp_minNeq=NaN;
Tpp_Jac=[];
Tpp_EIGEN=[];
Tpp_DOM_EIGEN=[];
Tpp_feasibleCrit=NaN;
Tpp_stableCrit=NaN;

Tpp_found=0;
Tpp_stepNR_Precision=NaN;

Tpp_HOPF=0;
Tpp_L0=NaN;
Tpp_L0_MATCONT=NaN;

Tpp_HOPF_NR_ASS=NaN;
Tpp_HOPF_MIN_SpecNR=NaN;

Tpp_TRSCRT=0;

Tpp_M=NaN;

%%%%%%%%%%%%%%%%%%%%%%%
%%%% STARTING step %%%%
%%%%%%%%%%%%%%%%%%%%%%%

%% growth
RT_step=RT_null+RT_change.*(stepNR-1);

%% interaction strengths
A_step=A_null+A_change.*(stepNR-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ANALYSE NON-TRIVIAL EQ %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% check stability
[Neq_step,Jac_step,EIGEN_step,DOM_EIGEN_step,feasibleCrit_step,stableCrit_step]=code.stability_functions.func_eigenLotVolt(NRspec,RT_step,A_step);

%% Tipping point info
if (feasibleCrit_step==0 || stableCrit_step==0) && Tpp_found==0
%     Neq_step
%     EIGEN_step
%     stepNR
    disp('initial step not feasible/unstable')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% APPLY CHANGE PRECISION %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RT_change_Precision=RT_change./NRsteps_fullPrecision;
A_change_Precision=A_change./NRsteps_fullPrecision;

stepNR_Precision=0;
while true
    
    %% stepNRs
    stepNR_Precision=stepNR_Precision+1;
    if stepNR_Precision>=NRsteps_fullPrecision+1
        break
    end

    %% growth
    RT_step_Precision=RT_step+RT_change_Precision.*(stepNR_Precision-1);

    %% interaction strengths
    A_step_Precision=A_step+A_change_Precision.*(stepNR_Precision-1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% ANALYSE NON-TRIVIAL EQ %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% check stability
    [Neq_step,Jac_step,EIGEN_step,DOM_EIGEN_step,feasibleCrit_step,stableCrit_step]=code.stability_functions.func_eigenLotVoltFeasTRSH0(NRspec,RT_step_Precision,A_step_Precision);
    
    %% Tipping point info
    if (feasibleCrit_step==0 || stableCrit_step==0) && Tpp_found==0
        
        Tpp_found=1;
        Tpp_stepNR_Precision=stepNR_Precision;
        
        %% eigenLotVolt info
        Tpp_Neq=Neq_step;
        Tpp_minNeq=min(Neq_step);
        Tpp_Jac=Jac_step;
        Tpp_EIGEN=EIGEN_step;
        Tpp_DOM_EIGEN=DOM_EIGEN_step;
        Tpp_feasibleCrit=feasibleCrit_step;
        Tpp_stableCrit=stableCrit_step;
        Tpp_M=((stepNR-1)./NRsteps_full)+(1/NRsteps_full).*((stepNR_Precision-1)./NRsteps_fullPrecision);
        
        %% Hopf bifurcation
        if length(Tpp_DOM_EIGEN)==2 && isreal(Tpp_DOM_EIGEN)==0 && Tpp_feasibleCrit==1
            
            Tpp_HOPF=1;
            
            %% Lyapunov coefficient
            if sum(Tpp_EIGEN>0)==2 %% only works when two eigenvalues>0
                [Tpp_L0, Tpp_L0_MATCONT]=code.stability_functions.func_Lyapunov_LV(A_step_Precision,Jac_step,NRspec);
            end
            
            %% Test NR of alternative stable statses
            if changeNR<=NRchangeASSAnalysis
                [Tpp_HOPF_NR_ASS,~,~,~,Tpp_HOPF_MIN_SpecNR,~,~] = code.stability_functions.func_analyse_ASSUBS_LV(resultFolder,NRspec,RT_step_Precision,A_step_Precision);
            end
            
        end
        
        %% Loss of feasiblity - check if Transcritical/EQ remaining spec is stable
        NRinv=NaN;
        feasibleCrit_TRSCRTSUB=0;
        stableCrit_TRSCRTSUB=0;
        if Tpp_feasibleCrit==0 && sum(Tpp_EIGEN>0)==1 && sum(Tpp_Neq<=0)==1

            %% get abundances of subset/check feasible
            SpecNRs_TRSCRTSUB=find(Tpp_Neq>0);
            A_TRSCRTSUB=A_step_Precision(SpecNRs_TRSCRTSUB,SpecNRs_TRSCRTSUB);
            RT_TRSCRTSUB=RT_step_Precision(SpecNRs_TRSCRTSUB,:);
            [Neq_TRSCRTSUB,~,~,~,feasibleCrit_TRSCRTSUB,stableCrit_TRSCRTSUB]=code.stability_functions.func_eigenLotVoltFeasTRSH0(NRspec-1,RT_TRSCRTSUB,A_TRSCRTSUB);
            
            %% stability in context of full community
            Neq_FULL_TRSCRTSUB=zeros(NRspec,1);
            Neq_FULL_TRSCRTSUB(SpecNRs_TRSCRTSUB,:)=Neq_TRSCRTSUB;
            
            %% potential invadors
            GROWTH=(Neq_FULL_TRSCRTSUB==0).*(RT_step_Precision+A_step_Precision*Neq_FULL_TRSCRTSUB);
            NRinv=sum(GROWTH>0);

            %% check stability of subset
            if feasibleCrit_TRSCRTSUB==1 && stableCrit_TRSCRTSUB==1 && NRinv==0
                Tpp_TRSCRT=1;
            end
            
        end
        
        %% display info
        %fprintf('%d - %d - Tipping point found at M=%.2f\n',NETnr,changeNR,(stepNR_Precision./NRsteps_full))
        
        if Tpp_feasibleCrit==0
            if Tpp_TRSCRT==1
                fprintf('%d - %d - Transcritical bicurcation at M=%.4f\n',NETnr,changeNR,Tpp_M)
            else
                fprintf('%d - %d - Loss of feasibility at M=%.4f\n',NETnr,changeNR,Tpp_M)
                %Tpp_Neq
                %Tpp_EIGEN
                %stableCrit_TRSCRTSUB
                %NRinv
                %pause
            end
        else
            if Tpp_HOPF==0
                fprintf('>>> %d - %d - STRANGE OUTPUT - FEASIBLE - NO HOPF!!!! M=%.4f\n',NETnr,changeNR,Tpp_M)
            else
                if Tpp_L0<0
                    fprintf('%d - %d - Supercritical HOPF (limit cycle), L0: %.4f at M=%.4f, NRASS=%d, minSpecNR=%d\n',NETnr,changeNR,Tpp_L0,Tpp_M,Tpp_HOPF_NR_ASS,Tpp_HOPF_MIN_SpecNR)
                elseif Tpp_L0>0
                    fprintf('%d - %d - Subcritical HOPF (instability), L0: %.4f at M=%.4f, NRASS=%d, minSpecNR=%d\n',NETnr,changeNR,Tpp_L0,Tpp_M,Tpp_HOPF_NR_ASS,Tpp_HOPF_MIN_SpecNR)
                else
                    fprintf('%d - %d - Unknown HOPF (likely sum(EIGEN>0)>2), L0: %.4f at M=%.4f, NRASS=%d, minSpecNR=%d\n',NETnr,changeNR,Tpp_L0,Tpp_M,Tpp_HOPF_NR_ASS,Tpp_HOPF_MIN_SpecNR)
                end
            end
        end
        
        %% change number NRsteps (shortens time series)
        NRsteps_fullPrecision=floor(Tpp_stepNR_Precision.*1.1);
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% SAVE INFO in DATASERIES %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    NTReq_series(stepNR_Precision,:)=Neq_step;
    NTR_EIGEN_series(stepNR_Precision,:)=EIGEN_step;
    NTR_DOM_EIGEN_series(stepNR_Precision,1)=DOM_EIGEN_step(1,1);
    
end

if Tpp_found==0
    fprintf('%d - %d - No Tipping point FOUND!!\n',NETnr,changeNR)
end