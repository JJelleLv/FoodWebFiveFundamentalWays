function [NR_ASS,EIGEN_found,Neq_found,COMB_found,MIN_SpecNR,MAXAb_Neq_found,SUBSIZE_NRinv_MAT] = func_analyse_ASSUBS_LV(resultFolder,NRspec,RT,A)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% FIND Alternative Stable Sub-networks %%%%%%%%%%%%%%%%
%%%%%% Multi-species variety of Yodzis and Innes, 1992 %%%%%
%%%%%% Code by Jelle Lever, 2015 %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% BOTH TYPE I and TYPE II are stable %%%%%%%%%%%%%%%%%%
%%%%%% LV analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% search for all ASS combinations %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:nearlySingularMatrix');

NR_ASS=0;

EIGEN_found={}; %% I don't know the size on beforehand...
Neq_found={};
COMB_found={};

MIN_SpecNR=0;
MAXAb_Neq_found=0;

SUBSIZE_NRinv_MAT=zeros(NRspec,NRspec+1);

%% go through all combinations and caluclate stability
for NRspec_eq_rev=1:NRspec
    NRspec_eq=NRspec+1-NRspec_eq_rev;
    
    allCombData=load(sprintf('%s%sALL_comb_S%d%sALL_comb_%d',resultFolder,filesep,NRspec,filesep,NRspec_eq));
    ALL_comb=allCombData.ALL_comb;
    NRcomb=allCombData.NRcomb;
    
    for CombNR=1:NRcomb
        
        %% create RT and A of the subset
        RT_test=RT(ALL_comb(CombNR,:),1);
        A_test=A(ALL_comb(CombNR,:),ALL_comb(CombNR,:));
        
        %% clear warnings
        lastwarn('')
        
        %% calculate abundances
        Neq_test=-(A_test^-1)*RT_test;
        
        %% if singular or badly scaled
        [warn_text, id_text]=lastwarn;
        SINGULAR=strcmp(id_text,'MATLAB:singularMatrix');
        NEARSINGULAR=strcmp(id_text,'MATLAB:nearlySingularMatrix');
        ILLCOND=strcmp(id_text,'MATLAB:illConditionedMatrix');
        
        MATRIX_COND=0;
        if SINGULAR==1 || NEARSINGULAR==1 || ILLCOND==1
            MATRIX_COND=1;
        end
        
        SUM_Neqtest=sum(Neq_test>=0);
        if SUM_Neqtest==NRspec_eq && MATRIX_COND==0
            
            %% stability of subset
            Jac_test=zeros(NRspec_eq,NRspec_eq);
            Jac_test(1:NRspec_eq+1:NRspec_eq*NRspec_eq)=RT_test+A_test*Neq_test;
            Jac_test=Jac_test+A_test.*(Neq_test*ones(1,NRspec_eq));
            
            if sum(sum(isnan(Jac_test)))==0 && sum(sum(isinf(Jac_test)))==0;
                
                EIGEN_test=eig(Jac_test);
                
                SUMeigen_test=sum(EIGEN_test<0); %% if all are smaller than 0, stable eq.
                
                if SUMeigen_test==NRspec_eq
                    
                    %% stability in context of full community
                    Neq=zeros(NRspec,1);
                    Neq(ALL_comb(CombNR,:),1)=Neq_test;
                    
                    %% potential invadors
                    GROWTH=(Neq==0).*(RT+A*Neq);
                    NRinv=sum(GROWTH>0);
                    
                    SUBSIZE_NRinv_MAT(NRspec_eq,NRinv+1)=SUBSIZE_NRinv_MAT(NRspec_eq,NRinv+1)+1;
                    
                    %% ASS if there are no invadors, NRinv==0
                    if NRinv==0
                        NR_ASS=NR_ASS+1;
                        
                        Jac=zeros(NRspec,NRspec);
                        Jac(1:NRspec+1:NRspec*NRspec)=RT+A*Neq;
                        Jac=Jac+A.*(Neq*ones(1,NRspec));
                        
                        EIGEN=eig(Jac);
                        
                        EIGEN_found{NR_ASS}=EIGEN;
                        Neq_found{NR_ASS}=Neq;
                        COMB_found{NR_ASS}=ALL_comb(CombNR,:);
                        MIN_SpecNR=NRspec_eq;
                        if max(Neq)>=MAXAb_Neq_found
                            MAXAb_Neq_found=max(Neq);
                        end
                    end
                end
            end
        end
    end
end

