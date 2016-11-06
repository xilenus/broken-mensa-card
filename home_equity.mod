//endogenous variables
var
c p_h h n_m R d n_u y z_h u_lo n_do n_s lambda_cc lambda_ns lambda_nm lambda_nu z_a g_hat G
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


fun = @(???F(epsilon_j)???) epsilon_j;
u_lo = OMEGA*g_hat(-1)*n_s(-1)*(PSI*???F(epsilon_do_bar)+integral(fun,epsilon_do_bar,inf)) + (1-OMEGA*g_hat(-1)*n_s(-1))*(PSI*???F(epsilon_bar)???+integral(fun,epsilon_bar,inf));

ALPHA*z_h/h - p_h/c + BETA*E*p_h(+1)/c(+1) + lambda_cc*n_m*CHI*p_h = 0;
1/c - lambda_cc + BETA*E(lambda_cc(+1)*(1-n_m(+1))-R/c(+1)) = 0;
BETA*(OMEGA*g_hat*(PSI*???F(epsilon_do_bar(+1)+integral(fun,epsilon_do_bar(+1),inf))-OMEGA*g_hat*(PSI*???F(epsilon_bar(+1))???+integral(fun,epsilon_bar(+1),inf))-lambda_nu*(1-g_hat-OMEGA*g_hat*(1-???F(epsilon_do_bar(+1)))-lambda_nm*OMEGA*g_hat*(???F(epsilon_do_bar(+1))-F(epsilon_j(+1))???)) = lambda_ns;
ZETA/c - lambda_cc*(CHI*p_h*h-d(-1)) - lambda_nm = 0;
PSI - epsilon_bar - lambda_nm = 0;
PSI - epsilon_do_bar + lambda_nu - lambda_nm = 0;
KAPPA - z_a*XI/c + lambda_ns*(1-RHO_u) + lambda_nu = 0;
ALPHA*z_h/h - p_h/c + BETA*E*p_h(+1)/c(+1) + lambda_cc*n_m*CHI*p_h = 0;
1/c - lambda_cc + BETA*E*(lambda_cc(+1)*(1-n_m(+1))-R/c(+1)) = 0;
ZETA/c - lambda_cc*(CHI*p_h*h-d(-1)) + epsilon_bar - PSI = 0;
KAPPA - z_a*XI/c + epsilon_do_bar - epsilon_bar + (1-RHO_u)*G = 0;
G = BETA*E*(-OMEGA*g_hat*integral(fun,epsilon_bar(+1),epsilon_do_bar(+1))+(epsilon_do_bar(+1)-epsilon_bar(+1))*(1-g_hat+OMEGA*g_hat*(1-???F(epsilon_do_bar(+1))???)+epsilon_bar(+1)*OMEGA*g_hat*(???F(epsilon_do_bar(+1))-F(epsilon_bar(+1))???)));
integral(fun,epsilon_bar(+1),epsilon_do_bar(+1)) = SIGMA/(2*pi)^(1/2)*(exp((-epsilon_bar(+1))^2/(2*SIGMA^2))-exp((-epsilon_do_bar(+1))^2/(2*SIGMA^2)));
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
