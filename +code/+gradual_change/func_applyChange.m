function [INT_list_t_REL_end,totFeedRate_end,RT_end,A_end,RT_change,A_change,Tpp_found,Tpp_stepNR,Tpp_Neq,Tpp_minNeq,Tpp_Jac,Tpp_EIGEN,Tpp_DOM_EIGEN, ...
Tpp_feasibleCrit,Tpp_stableCrit,Tpp_HOPF,Tpp_L0,Tpp_L0_MATCONT,NTReq_series,NTR_EIGEN_series,NTR_DOM_EIGEN_series,NRsteps]=func_applyChange(NETnr,changeNR,NRsteps_full,NRspec,TL_list,Spec_NRprey_count, ...
mr,mc,Ri,Ki,Tk,NRint_t,NRint_c,INT_list_c,INT_list_t_REL_null,totFeedRate_null,RT_null,A_null,ONLY_herb_carn,CHANGE_int,NEW_FRAC_int_nrnd,NEW_VAR_int,CHANGE_totFeedRate,NEW_RELtotFeedRate_min,NEW_RELtotFeedRate_max)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MAKE FINAL NETWORK %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% change relative interaction strengths
if CHANGE_int==1
    [INT_list_t_REL_end]=code.gradual_change.func_changeRelInt_t(NRspec,NRint_t,INT_list_t_REL_null,TL_list,NEW_FRAC_int_nrnd,NEW_VAR_int,ONLY_herb_carn);
else
    INT_list_t_REL_end=INT_list_t_REL_null;
end

%% change predator saturation
if CHANGE_totFeedRate==1
    [totFeedRate_end,INT_list_t_REL_end]=code.gradual_change.func_changeTotFeedRate(NRspec,NRint_t,INT_list_t_REL_end,TL_list,Spec_NRprey_count,totFeedRate_null,NEW_RELtotFeedRate_min,NEW_RELtotFeedRate_max,ONLY_herb_carn);
else
    totFeedRate_end=totFeedRate_null;
end

%% make final LotVolt network
[RT_end,A_end]=code.stability_functions.func_makeLotVolt(NRspec,mr,mc,Ri,Ki,totFeedRate_end,Tk,NRint_t,INT_list_t_REL_end,NRint_c,INT_list_c);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MAKE CHANGE MAT/VECT %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% change in growth
RT_change=(RT_end-RT_null)./NRsteps_full;

%% change in interaction strengths
A_change=(A_end-A_null)./NRsteps_full;

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
Tpp_stepNR=NaN;

Tpp_HOPF=0;
Tpp_L0=NaN;
Tpp_L0_MATCONT=NaN;

%%%%%%%%%%%%%%%%%%%%%%
%%%% APPLY CHANGE %%%%
%%%%%%%%%%%%%%%%%%%%%%

NRsteps=NRsteps_full;
stepNR=0;
while true
    
    %% stepNRs
    stepNR=stepNR+1;
    if stepNR>=NRsteps+1
        break
    end

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
        
        Tpp_found=1;
        Tpp_stepNR=stepNR;
        
        %% eigenLotVolt info
        Tpp_Neq=Neq_step;
        Tpp_minNeq=min(Neq_step);
        Tpp_Jac=Jac_step;
        Tpp_EIGEN=EIGEN_step;
        Tpp_DOM_EIGEN=DOM_EIGEN_step;
        Tpp_feasibleCrit=feasibleCrit_step;
        Tpp_stableCrit=stableCrit_step;
        
        %% Hopf bifurcation
        if length(Tpp_DOM_EIGEN)==2 && isreal(Tpp_DOM_EIGEN)==0 && Tpp_feasibleCrit==1
            
            Tpp_HOPF=1;
            if sum(Tpp_EIGEN>0)==2 %% only works when two eigenvalues>0
                [Tpp_L0, Tpp_L0_MATCONT]=code.stability_functions.func_Lyapunov_LV(A_step,Jac_step,NRspec);
            end
        end
        
%         %% display info
%         %fprintf('%d - %d - Tipping point found at M=%.2f\n',NETnr,changeNR,(stepNR./NRsteps_full))
%         if Tpp_feasibleCrit==0
%             fprintf('%d - %d - Loss of feasibility at M=%.2f\n',NETnr,changeNR,(stepNR./NRsteps_full))
%         else
%             if Tpp_HOPF==0
%                 fprintf('>>> %d - %d - STRANGE OUTPUT - FEASIBLE - NO HOPF!!!! M=%.2f, Neq_min=%.2f\n',NETnr,changeNR,((stepNR-1)./NRsteps_full),min(Tpp_Neq))
%             else
%                 if Tpp_L0<0
%                     fprintf('%d - %d - Supercritical HOPF (limit cycle), L0: %.4f at M=%.2f, Neq_min=%.2f\n',NETnr,changeNR,Tpp_L0,((stepNR-1)./NRsteps_full),min(Tpp_Neq))
%                 elseif Tpp_L0>0
%                     fprintf('%d - %d - Subcritical HOPF (instability), L0: %.4f at M=%.2f, Neq_min=%.2f\n',NETnr,changeNR,Tpp_L0,((stepNR-1)./NRsteps_full),min(Tpp_Neq))
%                 else
%                     fprintf('%d - %d - Unknown HOPF (likely sum(EIGEN>0)>2), L0: %.4f at M=%.2f, Neq_min=%.2f\n',NETnr,changeNR,Tpp_L0,((stepNR-1)./NRsteps_full),min(Tpp_Neq))
%                 end
%             end
%         end
        
        %% change number NRsteps (shortens time series)
        NRsteps=floor(Tpp_stepNR.*1.1);

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% SAVE INFO in DATASERIES %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    NTReq_series(stepNR,:)=Neq_step;
    NTR_EIGEN_series(stepNR,:)=EIGEN_step;
    NTR_DOM_EIGEN_series(stepNR,1)=DOM_EIGEN_step(1,1);
    
end

if Tpp_found==0
    fprintf('%d - %d - No Tipping point FOUND!!\n',NETnr,changeNR)
end