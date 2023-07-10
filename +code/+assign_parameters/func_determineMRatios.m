function [mRatioList_herb,mRatioList_carn,minMRatioList_herb,maxMRatioList_herb,meanMRatioList_herb,minMRatioList_carn,maxMRatioList_carn,meanMRatioList_carn]=func_determineMRatios(mr,mc,TL_list,NRint_uns,CIRCLE_Int_list,INT_list_NET_uns)

%% empty data
mRatioListNR_herb=0;
mRatioList_herb=[];
mRatioListNR_carn=0;
mRatioList_carn=[];

%% make list of body-mass ratios
M=mr+mc;
for IntNR_uns=1:NRint_uns
    if sum(CIRCLE_Int_list==IntNR_uns)==0 %% not on a circle
        PREY=INT_list_NET_uns(IntNR_uns,1);
        PRED=INT_list_NET_uns(IntNR_uns,2);
        mRatio=M(PRED,1)./M(PREY,1);
        if PREY~=PRED && (TL_list(PREY,1)==1) %% carnivore interaction
            mRatioListNR_herb=mRatioListNR_herb+1;
            mRatioList_herb(mRatioListNR_herb,1)=mRatio;
        end
        if PREY~=PRED && (TL_list(PREY,1)>1) %% carnivore interaction
            mRatioListNR_carn=mRatioListNR_carn+1;
            mRatioList_carn(mRatioListNR_carn,1)=mRatio;
        end
    end
end

%% gather info
minMRatioList_herb=min(mRatioList_herb);
maxMRatioList_herb=max(mRatioList_herb);
meanMRatioList_herb=mean(mRatioList_herb);
minMRatioList_carn=min(mRatioList_carn);
maxMRatioList_carn=max(mRatioList_carn);
meanMRatioList_carn=mean(mRatioList_carn);