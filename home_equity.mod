close all

//endogenous variables
var

//basic variables (15/15)
c h d
epsilon_bar epsilon_do_bar
n_m n_u n_s
m_hat   //number of meetings
g_hat   //probability of unemployed member get a job offer
g_f_hat //probability of firm fills a vacancy
v_hat   //number of vacancies
V       //firm value
p_h R

//additional variables (3/18)
G
F_bar F_do_bar
z_a z_h
;

//exogenous variables
varexo

//VAR(1) in logs (2)
eps_z_h //housing preference shock with mean 1
eps_z_a //productivity shock with mean 1
;

//parameters
parameters
pi GAMMA ALPHA_pat PSI SIGMA_eps ZETA RHO_u ETA MU THETA OMEGA XI KAPPA CHI PHI_h PHI_a Z_A Z_H
;

//initialization of parameters
pi = 3.1416;        //pi value
GAMMA = 0.9975;     //discount factor
ALPHA_pat = 0.043;  //housing preference
PSI = -7.068;       //utility from new location
SIGMA_eps = 3.47;   //standard deviation location preference shock
SIGMA_z_h = 0.26;   //standard deviation of housing preference shock
SIGMA_z_a = 0.013;  //standard deviation of productivity shock
ZETA = 1.6;         //moving cost
RHO_u = 0.035;      //rate of job destruction
ETA = 0.6;          //elasticity matching function
MU = 0.545;         //level parameter matching function
THETA = 0.181;      //vacancy cost
OMEGA = 0.3333;     //fraction of long-distance job offers
XI = 0.98;          //wage rule parameter
KAPPA = -4.381;     //utility from unemployment
CHI = 0.8;          //collateral requirement
PHI_h = 0.983;      //autocorrelation housing preference process
PHI_a = 0.983;      //autocorrelation technology process
Z_A = 0.9;          //shock coefficient productivity
Z_H = 0.9;          //shock coefficient housing preference


//model
model;

//conditions for the patient households (4/4)
c + p_h*(h-h(-1)) + ZETA*n_m +  R(-1)*d(-1) = (1 - n_u)*exp(z_a) - THETA*v_hat + d;
n_m = OMEGA*g_hat(-1)*n_s(-1)*F_do_bar + (1 - OMEGA*g_hat(-1)*n_s(-1))*F_bar;
n_s = n_u + RHO_u*(1 - n_u);
n_u = n_s(-1)*(1 - g_hat(-1) + OMEGA*g_hat(-1)*(1 - F_do_bar));

//first order conditions (4/8)
ALPHA_pat*exp(z_h)/h - p_h/c + GAMMA*p_h(+1)/c(+1) = 0;
1/c - GAMMA*R/c(+1) = 0;
ZETA/c + epsilon_bar - PSI = 0;
KAPPA + epsilon_do_bar - epsilon_bar + (1 - RHO_u)*G - exp(z_a)/c = 0;

//matching function (1/9)
m_hat = MU*n_s^ETA*v_hat^(1 - ETA);

//definitions for labor market aggregates (2/11)
g_hat = m_hat/n_s;
g_f_hat = m_hat/v_hat;

//equations for the firms (2/13)
V = (1 - XI)*exp(z_a) + (1 - RHO_u)*GAMMA*c/c(+1)*V(+1);
THETA = g_f_hat*(1 - OMEGA + OMEGA*F_do_bar(+1))*GAMMA*c/c(+1)*V(+1);

//market clearing conditions for housing and bonds (2/15)
h = 1;
d = 0;

//additional equations (3/18)
F_bar = 1/2*(1 + erf(epsilon_bar/(SIGMA_eps*2^0.5)));
F_do_bar = 1/2*(1 + erf(epsilon_do_bar/(SIGMA_eps*2^0.5)));
G = GAMMA*(-OMEGA*g_hat*(SIGMA_eps/(2*pi)^0.5*(exp(-epsilon_bar(+1)^2/(2*SIGMA_eps^2)) - exp(-epsilon_do_bar(+1)^2/(2*SIGMA_eps^2)))) +
(epsilon_do_bar(+1) - epsilon_bar(+1))*(1 - g_hat + OMEGA*g_hat*(1 - F_do_bar(+1))) + epsilon_bar(+1)*OMEGA*g_hat*(F_do_bar(+1) - F_bar(+1))) +
1/c*THETA*(MU/g_f_hat)^(1/ETA);

//shock AR(1)
z_a = -eps_z_a + Z_A*z_a(-1);
z_h = -eps_z_h + Z_H*z_h(-1);
end;

//assign steady state values
initval;
h = 1;
d = 0;
R = 1/GAMMA;
V = (XI - 1)/((1 - RHO_u)*GAMMA - 1);
G = 3.2897;
c = 0.9339;
g_hat = 0.6586;
v_hat = 0.1002;
epsilon_bar = PSI - ZETA/c;
F_bar = 1/2*(1 + erf(epsilon_bar/(SIGMA_eps*2^0.5)));
epsilon_do_bar = 1/c - KAPPA + (1 - RHO_u)*G + epsilon_bar;
F_do_bar = 1/2*(1 + erf(epsilon_do_bar/(SIGMA_eps*2^0.5)));
g_f_hat = THETA/((1 - OMEGA + OMEGA*F_do_bar)*GAMMA*V);
m_hat = g_f_hat*v_hat;
n_s = (m_hat/(MU*v_hat^(1 - ETA)))^(1/ETA);
g_hat = m_hat/n_s;
n_u = n_s*(1 - g_hat + OMEGA*g_hat*(1 - F_do_bar));
n_m = OMEGA*g_hat*n_s*F_do_bar + (1 - OMEGA*g_hat*n_s)*F_bar;
end;

//calc. and check steady state
steady;

//check eigenvalues
check;

//standard deviation of shocks
shocks;
var eps_z_a = SIGMA_z_a^2;
var eps_z_h = SIGMA_z_h^2;
end;

//stochastic simulation
stoch_simul(order=1,irf=25);
