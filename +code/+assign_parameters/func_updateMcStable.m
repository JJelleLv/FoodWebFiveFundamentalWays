function [mc,stableFound]=func_updateMcStable(mc,NRspec,TL_list,NRBASAL,NONBASAL_SpecNRs,NRint_t,INT_list_t_REL,thetaNeq_mat,thetaNeq,kappa,Neq_target,fj,aj,at,fr,ar,b,rKratio,MIN_predSat,MAX_predSat,PARAM_MAX_EIG)

%% predSatFound
stableFound=0;

%% determine predSat
[predSat,predSat_nonbs,totFeedRate,Jk,Tk]=code.assign_parameters.func_determinePredSat(NRspec,mc,fj,aj,at,b,thetaNeq_mat,thetaNeq,NONBASAL_SpecNRs);

%% check deviation from range
checkPredSat_nonbs=sum((predSat_nonbs>MAX_predSat)+(predSat_nonbs<MIN_predSat));
if checkPredSat_nonbs~=0
    error('initial mc has incorrect predSat')
end

%% determine mr, Ri and Ki, and LotVolt
[mr,Ri,Ki,RT,A,Jac,MATRIX_COND,EIGEN,stableCrit]=code.assign_parameters.func_assignMrRiKiLotVolt(NRspec,TL_list,NRBASAL,NRint_t,INT_list_t_REL,kappa,totFeedRate,Tk,Neq_target,fr,ar,b,rKratio,PARAM_MAX_EIG);
DOM_EIGEN=max(real(EIGEN));

%% return if already good
if stableCrit==1
    stableFound=1;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% update mc to get all in range %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ATTNR=0;
while true
    
    ATTNR=ATTNR+1;
    if ATTNR>=1e5
        %disp('stable failure')
        break
    end
    
    %% store old mc
    mc_old=mc;
    DOM_EIGEN_old=DOM_EIGEN;
        
    %% update mc
    dCh=0.1;
    chVect=1-dCh+2.*rand(NRspec,1).*dCh;
    mc=mc.*chVect;
    
    %% determine predSat
    [predSat,predSat_nonbs,totFeedRate,Jk,Tk]=code.assign_parameters.func_determinePredSat(NRspec,mc,fj,aj,at,b,thetaNeq_mat,thetaNeq,NONBASAL_SpecNRs);
    
    %% check deviation from range
    checkPredSat_nonbs=sum((predSat_nonbs>MAX_predSat)+(predSat_nonbs<MIN_predSat));
    
    %% check predSat
    if checkPredSat_nonbs~=0
        
        %fprintf('%d - undo change - wrong predSat\n',ATTNR)
        mc=mc_old;
        DOM_EIGEN=DOM_EIGEN_old;
        
    else
        
        %% determine mr, Ri and Ki, and LotVolt
        [mr,Ri,Ki,RT,A,Jac,MATRIX_COND,EIGEN,stableCrit]=code.assign_parameters.func_assignMrRiKiLotVolt(NRspec,TL_list,NRBASAL,NRint_t,INT_list_t_REL,kappa,totFeedRate,Tk,Neq_target,fr,ar,b,rKratio,PARAM_MAX_EIG);
        DOM_EIGEN=max(real(EIGEN));
        
        if MATRIX_COND==1
            mc=mc_old;
            DOM_EIGEN=DOM_EIGEN_old;
            %fprintf('%d - change undone, EIGEN: %.6f\n',ATTNR,DOM_EIGEN)
            %disp('MATRIX_COND==1')
        else
            if DOM_EIGEN_old<DOM_EIGEN %% may not be a number when MATRIX_CONC==1
                mc=mc_old;
                DOM_EIGEN=DOM_EIGEN_old;
                %fprintf('%d - change undone, EIGEN: %.6f\n',ATTNR,DOM_EIGEN)
            else
                %fprintf('%d - change accepted, EIGEN: %.6f\n',ATTNR,DOM_EIGEN)
                ATTNR=0;
            end
        end
        
    end
    
    if stableCrit==1
        stableFound=1;
        return
    end
    
end
