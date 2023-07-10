function [predSat,predSat_nonbs,totFeedRate,Jk,Tk]=func_determinePredSat(NRspec,mc,fj,aj,at,b,thetaNeq_mat,thetaNeq,NONBASAL_SpecNRs)

%% Maximum feeding rate of predators
Jk=fj.*aj.*(mc.^b);
Jk(Jk==0)=NaN;

%% respiration predators
Tk=at.*(mc.^b);
Tk(Tk==0)=NaN;

%% determine total feeding rate (thetaNeq_mat*totFeedRate-Tk=0)
thetaNeq_mat_nonbs=thetaNeq_mat(NONBASAL_SpecNRs,NONBASAL_SpecNRs);
Tk_nonbs=Tk(NONBASAL_SpecNRs,:);
totFeedRate_nonbs=((thetaNeq_mat_nonbs^-1)*Tk_nonbs);
totFeedRate=nan(NRspec,1);
totFeedRate(NONBASAL_SpecNRs,:)=totFeedRate_nonbs;

%% determine predator saturation
predSat=(totFeedRate.*thetaNeq)./Jk;
predSat_nonbs=predSat(NONBASAL_SpecNRs,:);