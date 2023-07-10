function [mc,mRatioFound]=func_updateMRatio(mc,NRspec,TL_list,NRBASAL,NONBASAL_SpecNRs,NRint_t,CIRCLE_Int_list,INT_list_t_REL,thetaNeq_mat,thetaNeq,kappa,Neq_target,fj,aj,at,fr,ar,b,rKratio,MIN_predSat,MAX_predSat,PARAM_MAX_EIG,MIN_mRatio_carn,MAX_mRatio_carn)

%% predSatFound
mRatioFound=0;

%% determine predSat
[predSat,predSat_nonbs,totFeedRate,Jk,Tk]=code.assign_parameters.func_determinePredSat(NRspec,mc,fj,aj,at,b,thetaNeq_mat,thetaNeq,NONBASAL_SpecNRs);

%% check deviation from range predSat
checkPredSat_nonbs=sum((predSat_nonbs>MAX_predSat)+(predSat_nonbs<MIN_predSat));
if checkPredSat_nonbs~=0
    error('initial mc has incorrect predSat')
end

%% determine mr, Ri and Ki, and LotVolt
[mr,Ri,Ki,RT,A,Jac,MATRIX_COND,EIGEN,stableCrit]=code.assign_parameters.func_assignMrRiKiLotVolt(NRspec,TL_list,NRBASAL,NRint_t,INT_list_t_REL,kappa,totFeedRate,Tk,Neq_target,fr,ar,b,rKratio,PARAM_MAX_EIG);
DOM_EIGEN=max(real(EIGEN));

%% check stability
if stableCrit==0
    error('initial mc is unstable')
end

%% determine mRatios
[mRatioList_herb,mRatioList_carn,minMRatioList_herb,maxMRatioList_herb,meanMRatioList_herb,minMRatioList_carn,maxMRatioList_carn,meanMRatioList_carn]=code.assign_parameters.func_determineMRatios(mr,mc,TL_list,NRint_t,CIRCLE_Int_list,INT_list_t_REL);

%% mRatio criterium
dMRatio_carn=(mRatioList_carn>MAX_mRatio_carn).*(mRatioList_carn-MAX_mRatio_carn)+(mRatioList_carn<MIN_mRatio_carn).*(MIN_mRatio_carn-mRatioList_carn);
sum_dMRatio_carn=sum(dMRatio_carn);

% fprintf('INITIAL mRatio:\n %.2f || %.2f - %.2f - %.2f || %.2f - %.2f - %.2f\n', ...
%                     sum_dMRatio_carn,meanMRatioList_carn,minMRatioList_carn,maxMRatioList_carn,meanMRatioList_herb,minMRatioList_herb,maxMRatioList_herb)

%% return if already good
if sum_dMRatio_carn==0
    mRatioFound=1;
    return
end

%% return when carn ratio is really bad
if sum_dMRatio_carn>1e4 %|| sum_dMRatio_herb>1e5
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% update mc to get all in range %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ATTNR=0;
while true
    
    ATTNR=ATTNR+1;
    if ATTNR>=1e5%1e5
        %disp('mRatio failure')
        break
    end
    
    %% store old mc
    mc_old=mc;
    sum_dMRatio_carn_old=sum_dMRatio_carn;

    %% update mc
    dCh=0.1;
    chVect=1-dCh+2.*rand(NRspec,1).*dCh;
    mc=mc.*chVect;
    
    %% determine predSat
    [~,predSat_nonbs,totFeedRate,~,Tk]=code.assign_parameters.func_determinePredSat(NRspec,mc,fj,aj,at,b,thetaNeq_mat,thetaNeq,NONBASAL_SpecNRs);
    
    %% check deviation from range
    checkPredSat_nonbs=sum((predSat_nonbs>MAX_predSat)+(predSat_nonbs<MIN_predSat));
    
    %% check predSat
    if checkPredSat_nonbs~=0
        %fprintf('%d - change undone (predSat), mRatio:\n %.2f || %.2f - %.2f - %.2f || %.2f - %.2f - %.2f\n', ...
        %    ATTNR,sum_dMRatio_carn,meanMRatioList_carn,minMRatioList_carn,maxMRatioList_carn,meanMRatioList_herb,minMRatioList_herb,maxMRatioList_herb)
        mc=mc_old;
        sum_dMRatio_carn=sum_dMRatio_carn_old;
    else
        %% determine mr, Ri and Ki, and LotVolt
        [mr,Ri,Ki,RT,A,Jac,MATRIX_COND,EIGEN,stableCrit]=code.assign_parameters.func_assignMrRiKiLotVolt(NRspec,TL_list,NRBASAL,NRint_t,INT_list_t_REL,kappa,totFeedRate,Tk,Neq_target,fr,ar,b,rKratio,PARAM_MAX_EIG);
        DOM_EIGEN=max(real(EIGEN));
        
        %% check stability
        if stableCrit==0
            %fprintf('%d - change undone (unstable), mRatio:\n %.2f || %.2f - %.2f - %.2f || %.2f - %.2f - %.2f\n', ...
            %    ATTNR,sum_dMRatio_carn,meanMRatioList_carn,minMRatioList_carn,maxMRatioList_carn,meanMRatioList_herb,minMRatioList_herb,maxMRatioList_herb)
            mc=mc_old;
            sum_dMRatio_carn=sum_dMRatio_carn_old;
        else
            %% determine mRatios
            [mRatioList_herb,mRatioList_carn,minMRatioList_herb,maxMRatioList_herb,meanMRatioList_herb,minMRatioList_carn,maxMRatioList_carn,meanMRatioList_carn]=code.assign_parameters.func_determineMRatios(mr,mc,TL_list,NRint_t,CIRCLE_Int_list,INT_list_t_REL);
            
            %% mRatio criterium
            dMRatio_carn=(mRatioList_carn>MAX_mRatio_carn).*(mRatioList_carn-MAX_mRatio_carn)+(mRatioList_carn<MIN_mRatio_carn).*(MIN_mRatio_carn-mRatioList_carn);
            sum_dMRatio_carn=sum(dMRatio_carn);
            
            %% check closer to desired mRatio
            if sum_dMRatio_carn_old<sum_dMRatio_carn
                %fprintf('%d - change undone (mRatio), mRatio:\n %.2f || %.2f - %.2f - %.2f || %.2f - %.2f - %.2f\n', ...
                %    ATTNR,sum_dMRatio_carn,meanMRatioList_carn,minMRatioList_carn,maxMRatioList_carn,meanMRatioList_herb,minMRatioList_herb,maxMRatioList_herb)
                mc=mc_old;
                sum_dMRatio_carn=sum_dMRatio_carn_old;
            elseif (sum_dMRatio_carn_old-sum_dMRatio_carn)>=0.001
                %fprintf('%d - change accepted, mRatio:\n %.2f || %.2f - %.2f - %.2f || %.2f - %.2f - %.2f\n', ...
                %    ATTNR,sum_dMRatio_carn,meanMRatioList_carn,minMRatioList_carn,maxMRatioList_carn,meanMRatioList_herb,minMRatioList_herb,maxMRatioList_herb)
                ATTNR=0;
            end
        end
    end
    
    if mod(ATTNR,1e4)==0 && ATTNR~=0
        %fprintf('%d - mRatio: %.2f \n',ATTNR,sum_dMRatio_carn)
    end
    
    %% return when critium fulfilled
    if sum_dMRatio_carn==0
        mRatioFound=1;
        return
    end
    
end
