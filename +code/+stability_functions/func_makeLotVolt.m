function [RT,A]=func_makeLotVolt(NRspec,mr,mc,Ri,Ki,totFeedRate,Tk,NRint_t,INT_list_t_REL,NRint_c,INT_list_c)

%%%%%%%%%%%%%%%%%
%%%% make RT %%%%
%%%%%%%%%%%%%%%%%

Ri(isnan(Ri))=0;
Tk(isnan(Tk))=0;
RT=Ri-Tk;

%%%%%%%%%%%%%%%%
%%%% make A %%%%
%%%%%%%%%%%%%%%%

%% make empty A
A=zeros(NRspec,NRspec);

%% add competitive interactions
for IntNR_c=1:NRint_c
    Spec1=INT_list_c(IntNR_c,1);
    Spec2=INT_list_c(IntNR_c,2);
    A(Spec1,Spec2)=-(Ri(Spec1,1).*INT_list_c(IntNR_c,3))./Ki(Spec1,1);
    A(Spec2,Spec1)=-(Ri(Spec2,1).*INT_list_c(IntNR_c,3))./Ki(Spec2,1);
end

%% add trophic interactions
for IntNR_t=1:NRint_t
    PREY=INT_list_t_REL(IntNR_t,1);
    PRED=INT_list_t_REL(IntNR_t,2);
    A(PRED,PREY)=A(PRED,PREY)+totFeedRate(PRED,1).*INT_list_t_REL(IntNR_t,3);
    A(PREY,PRED)=A(PREY,PRED)-(totFeedRate(PRED,1).*INT_list_t_REL(IntNR_t,3))./((1-INT_list_t_REL(IntNR_t,4)).*INT_list_t_REL(IntNR_t,5));
end