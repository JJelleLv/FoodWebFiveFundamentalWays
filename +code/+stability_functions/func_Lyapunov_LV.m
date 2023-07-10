function [L0, L0_MATCONT]=func_Lyapunov_LV(A,Jac,NRspec)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% function that calculates the Lyapunov coefficient %%%%%%%%%
%%%% for a symple Lotka-Volterra system dN/dt=RT*N+(A*N).*N %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Eigen vectors and values of Jac and adjoined Jac
[EIGEN_VECT,DIAG_EIGEN]=eig(Jac);
[ADJ_EIGEN_VECT,ADJ_DIAG_EIGEN]=eig(Jac');

%% the number of the eigenvalue with real >-1e-5 (this should be 0);
NR_EIG0=find(sum((real(DIAG_EIGEN)>0).*(imag(DIAG_EIGEN)>0))); %% imaginary >0
ADJ_NR_EIG0=find(sum((real(ADJ_DIAG_EIGEN)>0).*(imag(ADJ_DIAG_EIGEN)<0))); %% imaginary <0

%% q
q=EIGEN_VECT(:,NR_EIG0); %% matlab already scaled it so that q'*q=1;

%% p
p_adj_unsc=ADJ_EIGEN_VECT(:,ADJ_NR_EIG0);
p_scalefact=p_adj_unsc'*q;
p=p_adj_unsc/(conj(p_scalefact));

%% imaginary part of eigenvalue
w=imag(DIAG_EIGEN(NR_EIG0,NR_EIG0));

%% B derivatives - in LodkaVolterra the second order derivatives are equal to interaction matrix A!
B_DERIV_MAT_2D=A;
B_DERIV_MAT_2D(1:NRspec+1:NRspec*NRspec)=2*B_DERIV_MAT_2D(1:NRspec+1:NRspec*NRspec);

%% Bqq_conj, nondiagonal: q(k,1)*conj(q(l,1))+conj(q(k,1))*q(l,1), diagonal: q(k,1)*conj(q(k,1))
qq_conj_MAT_2D=q*transpose(conj(q))+conj(q)*transpose(q);
qq_conj_MAT_2D(1:NRspec+1:NRspec*NRspec)=q.*(conj(q));
Bqq_conj=sum(B_DERIV_MAT_2D.*qq_conj_MAT_2D,2);

%% ABqq_conj
ABqq_conj=(Jac^-1)*Bqq_conj;

%% BqABqq_conj
qABqq_conj_MAT_2D=q*transpose(ABqq_conj)+ABqq_conj*transpose(q);
qABqq_conj_MAT_2D(1:NRspec+1:NRspec*NRspec)=q.*ABqq_conj;
BqABqq_conj=sum(B_DERIV_MAT_2D.*qABqq_conj_MAT_2D,2);

%% Bqq
qq_MAT_2D=q*transpose(q)+q*transpose(q);
qq_MAT_2D(1:NRspec+1:NRspec*NRspec)=q.*(q);
Bqq=sum(B_DERIV_MAT_2D.*qq_MAT_2D,2);

%% iwIminA_Bqq
iwIminA=2*i*w*(eye(NRspec))-Jac;
iwIminA_Bqq=(iwIminA^-1)*Bqq;

%% Bq_conj_iwIminA_Bqq
q_conj_iwIminA_Bqq_MAT_2D=conj(q)*transpose(iwIminA_Bqq)+iwIminA_Bqq*transpose(conj(q));
q_conj_iwIminA_Bqq_MAT_2D(1:NRspec+1:NRspec*NRspec)=conj(q).*(iwIminA_Bqq);
Bq_conj_iwIminA_Bqq=sum(B_DERIV_MAT_2D.*q_conj_iwIminA_Bqq_MAT_2D,2);

%% Lyapunov coefficient
L0=(1./(2*w)).*real(-2*(p'*BqABqq_conj)+p'*Bq_conj_iwIminA_Bqq);
L0_MATCONT=(1./(2)).*real(-2*(p'*BqABqq_conj)+p'*Bq_conj_iwIminA_Bqq);
