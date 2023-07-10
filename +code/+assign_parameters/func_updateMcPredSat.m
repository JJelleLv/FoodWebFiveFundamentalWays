function [mc,predSatFound]=func_updateMcPredSat(mc,NRspec,NONBASAL_SpecNRs,thetaNeq_mat,thetaNeq,fj,aj,at,b,MIN_predSat,MAX_predSat)

%% predSatFound
predSatFound=0;

%% determine predSat
[~,predSat_nonbs,~,~,~]=code.assign_parameters.func_determinePredSat(NRspec,mc,fj,aj,at,b,thetaNeq_mat,thetaNeq,NONBASAL_SpecNRs);

%% total deviation from range
dPredSat_nonbs=(predSat_nonbs>MAX_predSat).*(predSat_nonbs-MAX_predSat)+(predSat_nonbs<MIN_predSat).*(MIN_predSat-predSat_nonbs);
sum_dPredSat_nonbs=sum(dPredSat_nonbs);

%% return if already good
if sum_dPredSat_nonbs==0
    predSatFound=1;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% update mc to get all in range %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ATTNR=0;
while true
    
    ATTNR=ATTNR+1;
    if ATTNR>=1e5
        %disp('predSat failure')
        break
    end
    
    %% store old mc
    mc_old=mc;
    sum_dPredSat_nonbs_old=sum_dPredSat_nonbs;
        
    %% update mc
    dCh=0.1;
    chVect=1-dCh+2.*rand(NRspec,1).*dCh;
    mc=mc.*chVect;
    
    %% determine predSat
    [predSat,predSat_nonbs,~,~,~]=code.assign_parameters.func_determinePredSat(NRspec,mc,fj,aj,at,b,thetaNeq_mat,thetaNeq,NONBASAL_SpecNRs);
    
    %% total deviation from range
    dPredSat_nonbs=(predSat_nonbs>MAX_predSat).*(predSat_nonbs-MAX_predSat)+(predSat_nonbs<MIN_predSat).*(MIN_predSat-predSat_nonbs);
    sum_dPredSat_nonbs=sum(dPredSat_nonbs);
    
    %% undo mc change
    if sum_dPredSat_nonbs>sum_dPredSat_nonbs_old
        mc=mc_old;
        sum_dPredSat_nonbs=sum_dPredSat_nonbs_old;
        %fprintf('%d - change undone, sum_dPredSat: %.4f\n',ATTNR,sum_dPredSat_nonbs)
    elseif (sum_dPredSat_nonbs_old-sum_dPredSat_nonbs)>=0.001
        %% reset counter, decrease in sum_dPredSat_nonbs_old
        ATTNR=0;
        %fprintf('%d - change accepted, sum_dPredSat SMALLER: %.4f\n',ATTNR,sum_dPredSat_nonbs)
    else
        %fprintf('%d - change accepted, sum_dPredSat EQUAL: %.4f\n',ATTNR,sum_dPredSat_nonbs)
    end
    
    if sum_dPredSat_nonbs==0
        predSatFound=1;
        break
    end
end
