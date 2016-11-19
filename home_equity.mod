//endogenous variables
var
c p_h h n_m R d n_u y u_lo n_do n_s lambda_cc lambda_ns lambda_nm lambda_nu z_a G
c_squiggle h_squiggle d_squiggle n_m_squiggle n_u_squiggle n_s_squiggle n_do_squiggle
E epsilon_bar epsilon_do_bar
epsilon_bar_squiggle epsilon_do_bar_squiggle
g_hat n_s_hat n_u_hat m_hat g_f_hat v_hat V
lambda_ns_squiggle lambda_nm_squiggle lambda_nu_squiggle G_squiggle
lambda_t_squiggle
f_do_bar F_do_bar f_bar F_bar f_j F_j f_do_bar_squiggle F_do_bar_squiggle f_bar_squiggle F_bar_squiggle
;

//exogenous variables
varexo
epsilon_j z_h
;

//parameters
parameters
BETA GAMMA ALPHA_imp ALPHA_pat PSI SIGMA ZETA RHO_u ETA MU THETA OMEGA XI KAPPA NU CHI PHI_h PHI_a pi expon
;

//initialization of parameters
pi = 3.14159265359;
expon = 2.71828182846;
BETA = 0.9899;
GAMMA = 0.9975;
ALPHA_imp = 0.14;
ALPHA_pat = 0.043;
PSI = -7.068;
SIGMA = 3.47;
ZETA = 1.6;
RHO_u = 0.035;
ETA = 0.6;
MU = 0.545;
THETA = 0.181;
OMEGA = 1/3;
XI = 0.98;
KAPPA = -4.381;
NU = 0.2;
CHI = 0.8;
PHI_h = 0.983;
PHI_a = 0.983;

//model
model;
//INTEGRAL EQUATION
//epsilon^2*(1+erf(epsilon/(SIGMA*2^0.5)))/4(HIGHER BOUND) - epsilon^2*(1+erf(epsilon/(SIGMA*2^0.5)))/4(LOWER BOUND)

//optimization problem for impatient households

f_do_bar = 1/(2*SIGMA^2*pi)^0.5*expon^(-epsilon_do_bar^2/(2*SIGMA^2));
F_do_bar = (1 + erf(epsilon_do_bar/(SIGMA*2^0.5)))/2;
f_bar = 1/(2*SIGMA^2*pi)^0.5*expon^(-epsilon_bar^2/(2*SIGMA^2));
F_bar = (1 + erf(epsilon_bar/(SIGMA*2^0.5)))/2;
f_j = 1/(2*SIGMA^2*pi)^0.5*expon^(-epsilon_j^2/(2*SIGMA^2));
F_j = (1 + erf(epsilon_j/(SIGMA*2^0.5)))/2;

c + p_h*(h-h(-1)) + ZETA*n_m + R(-1)*d(-1) = (1-n_u)*z_a*XI + d;
d = n_m*CHI*p_h*h + (1-n_m)*d(-1);
n_s = n_u + RHO_u*(1-n_u);
n_u = n_s(-1)*(1-g_hat(-1)+OMEGA*g_hat(-1)*(1-F_do_bar));
n_m = OMEGA*g_hat(-1)*n_s(-1)*F_do_bar + (1-OMEGA*g_hat(-1)*n_s(-1))*F_bar;




u_lo = OMEGA*g_hat(-1)*n_s(-1)*(PSI*F_do_bar+(1 - epsilon_do_bar^2*(1+erf(epsilon_do_bar/(SIGMA*2^0.5)))/4)) + (1-OMEGA*g_hat(-1)*n_s(-1))*(PSI*F_bar+(1 - epsilon_bar^2*(1+erf(epsilon_bar/(SIGMA*2^0.5)))/4));

//first order conditions for impatient households

ALPHA_imp*z_h/h - p_h/c + BETA*E*p_h(+1)/c(+1) + lambda_cc*n_m*CHI*p_h = 0;
1/c - lambda_cc + BETA*E*(lambda_cc(+1)*(1-n_m(+1))-R/c(+1)) = 0;

lambda_ns = BETA*(OMEGA*g_hat*(PSI*F_do_bar(+1)+(1 - epsilon_do_bar(+1)^2*(1+erf(epsilon_do_bar(+1)/(SIGMA*2^0.5)))/4))-OMEGA*g_hat*(PSI*F_bar(+1)+(1 - epsilon_bar(+1)^2*(1+erf(epsilon_bar(+1)/(SIGMA*2^0.5)))/4))-lambda_nu*(1-g_hat-OMEGA*g_hat*(1-F_do_bar(+1))-lambda_nm*OMEGA*g_hat*(F_do_bar(+1)-F_j(+1))));

ZETA/c - lambda_cc*(CHI*p_h*h-d(-1)) - lambda_nm = 0;
PSI - epsilon_bar - lambda_nm = 0;
PSI - epsilon_do_bar + lambda_nu - lambda_nm = 0;
KAPPA - z_a*XI/c + lambda_ns*(1-RHO_u) + lambda_nu = 0;

ALPHA_imp*z_h/h - p_h/c + BETA*E*p_h(+1)/c(+1) + lambda_cc*n_m*CHI*p_h = 0;
1/c - lambda_cc + BETA*E*(lambda_cc(+1)*(1-n_m(+1))-R/c(+1)) = 0;
ZETA/c - lambda_cc*(CHI*p_h*h-d(-1)) + epsilon_bar - PSI = 0;
KAPPA - z_a*XI/c + epsilon_do_bar - epsilon_bar + (1-RHO_u)*G = 0;

G = BETA*E*(-OMEGA*g_hat*(epsilon_do_bar(+1)^2*(1+erf(epsilon_do_bar(+1)/(SIGMA*2^0.5)))/4 - epsilon_bar(+1)^2*(1+erf(epsilon_bar(+1)/(SIGMA*2^0.5)))/4)+(epsilon_do_bar(+1)-epsilon_bar(+1))*(1-g_hat+OMEGA*g_hat*(1-F_do_bar(+1)))+epsilon_bar(+1)*OMEGA*g_hat*(F_do_bar(+1)-F_bar(+1)));

(epsilon_do_bar(+1)^2*(1+erf(epsilon_do_bar(+1)/(SIGMA*2^0.5)))/4 - epsilon_bar(+1)^2*(1+erf(epsilon_bar(+1)/(SIGMA*2^0.5)))/4) = SIGMA/(2*pi)^(1/2)*(exp((-epsilon_bar(+1))^2/(2*SIGMA^2))-exp((-epsilon_do_bar(+1))^2/(2*SIGMA^2)));


//optimization problem for the patient households

f_do_bar_squiggle = 1/(2*SIGMA^2*pi)^0.5*exp(-epsilon_do_bar_squiggle^2/(2*SIGMA^2));
F_do_bar_squiggle = (1 + erf(epsilon_do_bar_squiggle/(SIGMA*2^0.5)))/2;
f_bar_squiggle = 1/(2*SIGMA^2*pi)^0.5*exp(-epsilon_bar_squiggle^2/(2*SIGMA^2));
F_bar_squiggle = (1 + erf(epsilon_bar_squiggle/(SIGMA*2^0.5)))/2;

c_squiggle + p_h*(h_squiggle-h_squiggle(-1)) + ZETA*n_m_squiggle + R(-1)*d_squiggle(-1) = (1-n_u_squiggle)*z_a*XI + d_squiggle;
d_squiggle = n_m_squiggle*CHI*p_h*h_squiggle + (1-n_m_squiggle)*d_squiggle(-1);
n_s_squiggle = n_u_squiggle + RHO_u*(1-n_u_squiggle);
n_u_squiggle = n_s_squiggle(-1)*(1-g_hat(-1)+OMEGA*g_hat(-1)*(1-F_do_bar_squiggle));
n_m_squiggle = OMEGA*g_hat(-1)*n_s_squiggle(-1)*F_do_bar_squiggle + (1-OMEGA*g_hat(-1)*n_s_squiggle(-1))*F_bar_squiggle;




u_lo = OMEGA*g_hat(-1)*n_s_squiggle(-1)*(PSI*F_do_bar_squiggle+(1 - epsilon_do_bar_squiggle^2*(1+erf(epsilon_do_bar_squiggle/(SIGMA*2^0.5)))/4)) + (1-OMEGA*g_hat(-1)*n_s_squiggle(-1))*(PSI*F_bar_squiggle+(1 - epsilon_bar_squiggle^2*(1+erf(epsilon_bar_squiggle/(SIGMA*2^0.5)))/4));

//first order conditions for patient households

ALPHA_pat*z_h/h_squiggle - p_h/c_squiggle + GAMMA*E*p_h(+1)/c_squiggle(+1) = 0;
1/c_squiggle + GAMMA*E*(-R/c_squiggle(+1)) = 0;

GAMMA*(OMEGA*g_hat*(PSI*F_do_bar_squiggle(+1)+(1 - epsilon_do_bar_squiggle(+1)^2*(1+erf(epsilon_do_bar_squiggle(+1)/(SIGMA*2^0.5)))/4))-OMEGA*g_hat*(PSI*F_bar_squiggle(+1)+(1 - epsilon_bar_squiggle(+1)^2*(1+erf(epsilon_bar_squiggle(+1)/(SIGMA*2^0.5)))/4))-lambda_nu_squiggle*(1-g_hat-OMEGA*g_hat*(1-F_do_bar_squiggle(+1))-lambda_nm_squiggle*OMEGA*g_hat*(F_do_bar_squiggle(+1)-F_j(+1)))) = lambda_ns_squiggle;

ZETA/c_squiggle - lambda_nm_squiggle = 0;
PSI - epsilon_bar_squiggle - lambda_nm_squiggle = 0;
PSI - epsilon_do_bar_squiggle + lambda_nu_squiggle - lambda_nm_squiggle = 0;
KAPPA - z_a*XI/c_squiggle + lambda_ns_squiggle*(1-RHO_u) + lambda_nu_squiggle = 0;

ALPHA_pat*z_h/h_squiggle - p_h/c_squiggle + GAMMA*E*p_h(+1)/c_squiggle(+1) = 0;
1/c_squiggle + GAMMA*E*(-R/c_squiggle(+1)) = 0;
ZETA/c_squiggle + epsilon_bar_squiggle - PSI = 0;
KAPPA - z_a*XI/c_squiggle + epsilon_do_bar_squiggle - epsilon_bar_squiggle + (1-RHO_u)*G_squiggle = 0;

G_squiggle = GAMMA*E*(-OMEGA*g_hat*(epsilon_do_bar_squiggle(+1)^2*(1+erf(epsilon_do_bar_squiggle(+1)/(SIGMA*2^0.5)))/4 - epsilon_bar_squiggle(+1)^2*(1+erf(epsilon_bar_squiggle(+1)/(SIGMA*2^0.5)))/4)+(epsilon_do_bar_squiggle(+1)-epsilon_bar_squiggle(+1))*(1-g_hat+OMEGA*g_hat*(1-F_do_bar_squiggle(+1)))+epsilon_bar_squiggle(+1)*OMEGA*g_hat*(F_do_bar_squiggle(+1)-F_bar_squiggle(+1)));

(epsilon_do_bar_squiggle(+1)^2*(1+erf(epsilon_do_bar_squiggle(+1)/(SIGMA*2^0.5)))/4 - epsilon_bar_squiggle(+1)^2*(1+erf(epsilon_bar_squiggle(+1)/(SIGMA*2^0.5)))/4) = SIGMA/(2*pi)^(1/2)*(exp((-epsilon_bar_squiggle(+1))^2/(2*SIGMA^2))-exp((-epsilon_do_bar_squiggle(+1))^2/(2*SIGMA^2)));


epsilon_bar_squiggle = PSI - ZETA/c_squiggle;



n_u_hat = NU*n_u + (1 - NU)*n_u_squiggle;
n_s_hat = NU*n_s + (1 - NU)*n_s_squiggle;
m_hat = MU*n_s_hat^ETA*v_hat^(1-ETA);
g_hat = m_hat/n_s_hat;
g_f_hat = m_hat/v_hat;

y = XI*z_a;
V = (1 - XI)*z_a + (1 - RHO_u)*E*lambda_t_squiggle*V(+1);
lambda_t_squiggle = GAMMA*c_squiggle/c_squiggle(+1);

THETA = g_f_hat*E*(1-OMEGA+OMEGA*NU*n_s/n_s_hat*F_do_bar(+1)+OMEGA*(1-NU)*n_s_squiggle/n_s_hat*F_do_bar_squiggle(+1))*lambda_t_squiggle*V(+1);

//flow equations
n_m = n_do*F_do_bar+(1-n_do)*F_bar;
n_do = OMEGA*g_hat(-1)*n_s;
n_s = n_u+RHO_u*(1-n_u);
n_u = n_s(-1)*(1-g_hat(-1)+OMEGA*g_hat(-1)*(1-F_do_bar));


end;


//assign steady state values
initval;

n_u = 0.05;
m_hat = 0.0065;
F_do_bar = 0.001 - F_do_bar_squiggle;
h = 1.4*y/p_h;
h = h_squiggle;
g_f_hat = 0.34;

end;

//calc. and check steady state
steady;

//check eigenvalues
check;

//standard deviation of shocks
shocks;
var epsilon_j = SIGMA^2;
var z_h = SIGMA^2;
end;

//stochastic simulation
stoch_simul(order=1,irf=40);
