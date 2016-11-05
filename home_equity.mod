# broken-mensa-card

//endogenous variables
var
c p_h h n_m R d n_u y z_h u_lo n_do n_s lambda_cc z_a g_hat
E epsilon_bar epsilon_do_bar
;

//exogenous variables
varexo
epsilon_j
;

//parameters
parameters
BETA GAMMA ALPHA_imp ALPHA_pat PSI SIGMA ZETA RHO_u ETA MU THETA OMEGA XI KAPPA NU CHI PHI_h PHI_a
;

//initialization of parameters
ALPHA = ALPHA_imp;


//model
model;
c + p_h*(h-h(-1)) + ZETA*n_m + R(-1)*d(-1) = (1-n_u)*z_a*XI + d;
d = n_m*CHI*p_h*h + (1-n_m)*d(-1);
n_s = n_u + RHO_u*(1-n_u);
n_u = n_s(-1)*(1-g_hat(-1)+OMEGA*g_hat(-1)*(1-???F(epsilon_do_bar)???));
n_m = OMEGA*g_hat(-1)*n_s(-1)*???F(epsilon_do_bar)??? + (1-OMEGA*g_hat(-1)*n_s(-1))*???F(epsilon_bar)???;

p_h/c = ALPHA*z_h/h + BETA*E*p_h(+1)/c(+1) + lambda_cc*n_m*CHI*p_h;
1/c = BETA*E*(R/c(+1)-lambda_cc(+1)*(1-n_m(+1))) + lambda_cc;
u_lo = n_do*(PSI*)
end;

//assign steady state values
initval;

end;

//calc. and check steady state
steady;

//check eigenvalues
check;

//standard deviation of shocks
shocks;

end;

//stochastic simulation
stoch_simul();
