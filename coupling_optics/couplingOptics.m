% Fibre coupling params

w0=10e-6;

% D/f optimum from coupling_Ruilier (depends on lambda, w0, and D)

D_f = 0.135;



% f= tel.D/D_f;
f = 0.7/D_f;
 

EF_atm = comm_signal.amplitude.*exp(1i.*comm_signal.meanRmPhase);

FCE =  eta_Chen(EF_atm,w0,lambda,f,tel.D);