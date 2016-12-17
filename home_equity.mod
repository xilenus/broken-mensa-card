//endogenous variables
var
c h d 
c_squiggle h_squiggle d_squiggle

epsilon_bar epsilon_do_bar
epsilon_bar_squiggle epsilon_do_bar_squiggle

n_m n_u n_s n_do
n_m_squiggle n_u_squiggle n_s_squiggle n_do_squiggle

n_u_hat n_s_hat m_hat g_hat g_f_hat v_hat

V
p_h R
lambda_cc

f_do_bar F_do_bar f_bar F_bar
f_do_bar_squiggle F_do_bar_squiggle f_bar_squiggle F_bar_squiggle
f_j F_j

lambda_t_squiggle y z_a

G G_squiggle
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
OMEGA = 0.3333;
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
ALPHA_imp*z_h/h - p_h/c + BETA*p_h(+1)/c(+1) + lambda_cc*n_m*CHI*p_h = 0;
1/c - lambda_cc + BETA*(lambda_cc(+1)*(1-n_m(+1))-R/c(+1)) = 0;
PSI = ZETA/c + epsilon_bar - lambda_cc*(CHI*p_h*h-d(-1));
epsilon_do_bar - epsilon_bar = y/c - KAPPA + (1-RHO_u)*G;

n_s = n_u + RHO_u*(1-n_u);
n_u = n_s(-1)*(1-g_hat(-1)+OMEGA*g_hat(-1)*(1-F_do_bar));
n_m = OMEGA*g_hat(-1)*n_s(-1)*F_do_bar + (1-OMEGA*g_hat(-1)*n_s(-1))*F_bar;
n_do = OMEGA*g_hat(-1)*n_s;

//optimization problem for the patient households

f_do_bar_squiggle = 1/(2*SIGMA^2*pi)^0.5*exp(-epsilon_do_bar_squiggle^2/(2*SIGMA^2));
F_do_bar_squiggle = (1 + erf(epsilon_do_bar_squiggle/(SIGMA*2^0.5)))/2;
f_bar_squiggle = 1/(2*SIGMA^2*pi)^0.5*exp(-epsilon_bar_squiggle^2/(2*SIGMA^2));
F_bar_squiggle = (1 + erf(epsilon_bar_squiggle/(SIGMA*2^0.5)))/2;

c_squiggle + p_h*(h_squiggle-h_squiggle(-1)) + ZETA*n_m_squiggle + R(-1)*d_squiggle(-1) = (1-n_u_squiggle)*z_a*XI + d_squiggle;
d_squiggle = n_m_squiggle*CHI*p_h*h_squiggle + (1-n_m_squiggle)*d_squiggle(-1);
ALPHA_pat*z_h/h_squiggle - p_h/c_squiggle + GAMMA*p_h(+1)/c_squiggle(+1) = 0;
1/c_squiggle + GAMMA*(-R/c_squiggle(+1)) = 0;
PSI = ZETA/c_squiggle + epsilon_bar_squiggle - lambda_cc*(CHI*p_h*h_squiggle-d_squiggle(-1));
KAPPA - z_a*XI/c_squiggle + epsilon_do_bar_squiggle - epsilon_bar_squiggle + (1-RHO_u)*G_squiggle = 0;

n_s_squiggle = n_u_squiggle + RHO_u*(1-n_u_squiggle);
n_u_squiggle = n_s_squiggle(-1)*(1-g_hat(-1)+OMEGA*g_hat(-1)*(1-F_do_bar_squiggle));
n_m_squiggle = OMEGA*g_hat(-1)*n_s_squiggle(-1)*F_do_bar_squiggle + (1-OMEGA*g_hat(-1)*n_s_squiggle(-1))*F_bar_squiggle;
n_do_squiggle = OMEGA*g_hat(-1)*n_s_squiggle;

//LABOR MARKET

n_u_hat = NU*n_u + (1 - NU)*n_u_squiggle;
n_s_hat = NU*n_s + (1 - NU)*n_s_squiggle;
m_hat = MU*n_s_hat^ETA*v_hat^(1-ETA);
g_hat = m_hat/n_s_hat;
g_f_hat = m_hat/v_hat;

//FIRMS

V = (1 - XI)*z_a + (1 - RHO_u)*lambda_t_squiggle*V(+1);
THETA = g_f_hat*(1-OMEGA+OMEGA*NU*n_s/n_s_hat*F_do_bar(+1)+OMEGA*(1-NU)*n_s_squiggle/n_s_hat*F_do_bar_squiggle(+1))*lambda_t_squiggle*V(+1);

NU*h + (1-NU)*h_squiggle = 1;
NU*d + (1-NU)*d_squiggle = 0;

//OTHER EQUATIONS

lambda_t_squiggle = GAMMA*c_squiggle/c_squiggle(+1);
y = XI*z_a;
G = BETA*(-OMEGA*g_hat*(epsilon_do_bar(+1)^2*(1+erf(epsilon_do_bar(+1)/(SIGMA*2^0.5)))/4 - epsilon_bar(+1)^2*(1+erf(epsilon_bar(+1)/(SIGMA*2^0.5)))/4)+(epsilon_do_bar(+1)-epsilon_bar(+1))*(1-g_hat+OMEGA*g_hat*(1-F_do_bar(+1)))+epsilon_bar(+1)*OMEGA*g_hat*(F_do_bar(+1)-F_bar(+1)));
G_squiggle = GAMMA*(-OMEGA*g_hat*(epsilon_do_bar_squiggle(+1)^2*(1+erf(epsilon_do_bar_squiggle(+1)/(SIGMA*2^0.5)))/4 - epsilon_bar_squiggle(+1)^2*(1+erf(epsilon_bar_squiggle(+1)/(SIGMA*2^0.5)))/4)+(epsilon_do_bar_squiggle(+1)-epsilon_bar_squiggle(+1))*(1-g_hat+OMEGA*g_hat*(1-F_do_bar_squiggle(+1)))+epsilon_bar_squiggle(+1)*OMEGA*g_hat*(F_do_bar_squiggle(+1)-F_bar_squiggle(+1)));
//KAPPA - z_a*XI/c + epsilon_do_bar - epsilon_bar + (1-RHO_u)*G = 0;
//KAPPA - z_a*XI/c_squiggle + epsilon_do_bar_squiggle - epsilon_bar_squiggle + (1-RHO_u)*G_squiggle = 0;

end;









//assign steady state values
initval;






//GIVEN BY THE AUTHOR

n_u_hat = 0.05;          //First, the aggregate unemployment rate in the steady state is five percent
m_hat = 0.0065;          //Second, the steady-state aggregate mobility rate is 0.65 percent per month
g_f_hat = 0.34;          //Sixth, the probability that a vacancy is filled is 0.34 in the steady state


//WRITE YOUR OWN VALUES (THE ONLY INPUT VALUES)

h = 1; //stock of housing //fixed and normalized to one

R = 1.05; //gross interest rate on debt
g_hat = 0.25; //probability of unemployed member to get a job offer
n_s = 0.15; //fraction of impatient households that look for job
n_s_squiggle = 0.15; //fraction of patient households that look for job

z_a = 300; //worker produce per period
d = 2000; //amount of debt

c = 100; //non-durable consumption of impatient households
c_squiggle = 120; //non-durable consumption of patient households






// Calculate analytical steady state values of endog. variables, ratios etc



epsilon_j = 0;
f_j = 1/(2*SIGMA^2*pi)^0.5*expon^(-epsilon_j^2/(2*SIGMA^2));
F_j = (1 + erf(epsilon_j/(SIGMA*2^0.5)))/2;

h_squiggle = h;                              //Fifth, the credit-constrained households consume the same amount of housing in the steady state as the patient households
lambda_t_squiggle = GAMMA;
V = (1-XI)*z_a/(RHO_u*lambda_t_squiggle);
n_s_hat = n_s*NU + n_s_squiggle*(1-NU);
d_squiggle = NU*d/(NU-1);
v_hat = m_hat/g_f_hat;
y = XI*z_a;
n_do = OMEGA*g_hat*n_s;
p_h = 1.4*z_a*12/(h+h_squiggle);                            //Fourth, the steady-state value of housing wealth is 140 per cent of annual output
n_u = (n_s-RHO_u)/(1-RHO_u);
n_u_squiggle = (NU*n_u-n_u_hat)/(NU-1);


F_do_bar = (n_s-n_s*g_hat+OMEGA*n_s*g_hat-n_u)/(OMEGA*n_s*g_hat);
epsilon_do_bar = SIGMA*2^0.5*erfinv(2*F_do_bar - 1);
f_do_bar = 1/(2*SIGMA^2*pi)^0.5*expon^(-epsilon_do_bar^2/(2*SIGMA^2));


epsilon_do_bar_squiggle = erfinv(-(erf(epsilon_do_bar/(SIGMA*2^0.5))*NU+1.998)/(1-NU))*SIGMA*2^0.5;        //Third, the steady-state mobility rate due to members with long-distance job offers is 0.10 per cent per month
f_do_bar_squiggle = 1/(2*SIGMA^2*pi)^0.5*exp(-epsilon_do_bar_squiggle^2/(2*SIGMA^2));            //Third, the steady-state mobility rate due to members with long-distance job offers is 0.10 per cent per month
F_do_bar_squiggle = (1 + erf(epsilon_do_bar_squiggle/(SIGMA*2^0.5)))/2;                          //Third, the steady-state mobility rate due to members with long-distance job offers is 0.10 per cent per month

n_do_squiggle = OMEGA*g_hat*n_s_squiggle;

n_m = (z_a*XI-n_u*z_a*XI+d-c-R*d)/ZETA;
n_m_squiggle = (z_a*XI-n_u_squiggle*z_a*XI+d_squiggle-c_squiggle-R*d_squiggle)/ZETA;

F_bar = (n_m - OMEGA*g_hat*n_s*F_do_bar)/(1-OMEGA*g_hat*n_s);
epsilon_bar = SIGMA*2^0.5*erfinv(2*F_bar - 1);
f_bar = 1/(2*SIGMA^2*pi)^0.5*expon^(-epsilon_bar^2/(2*SIGMA^2));
F_bar_squiggle = (n_m_squiggle-OMEGA*g_hat*n_s_squiggle*F_do_bar_squiggle)/(1-OMEGA*g_hat*n_s_squiggle);
epsilon_bar_squiggle = SIGMA*2^0.5*erfinv(2*F_bar_squiggle - 1);
f_bar_squiggle = 1/(2*SIGMA^2*pi)^0.5*expon^(-epsilon_bar_squiggle^2/(2*SIGMA^2));
G = BETA*(-OMEGA*g_hat*(epsilon_do_bar^2*(1+erf(epsilon_do_bar/(SIGMA*2^0.5)))/4 - epsilon_bar^2*(1+erf(epsilon_bar/(SIGMA*2^0.5)))/4)+(epsilon_do_bar-epsilon_bar)*(1-g_hat+OMEGA*g_hat*(1-F_do_bar))+epsilon_bar*OMEGA*g_hat*(F_do_bar-F_bar));
G_squiggle = GAMMA*(-OMEGA*g_hat*(epsilon_do_bar_squiggle^2*(1+erf(epsilon_do_bar_squiggle/(SIGMA*2^0.5)))/4 - epsilon_bar_squiggle^2*(1+erf(epsilon_bar_squiggle/(SIGMA*2^0.5)))/4)+(epsilon_do_bar_squiggle-epsilon_bar_squiggle)*(1-g_hat+OMEGA*g_hat*(1-F_do_bar_squiggle))+epsilon_bar_squiggle*OMEGA*g_hat*(F_do_bar_squiggle-F_bar_squiggle));
lambda_cc = (1-BETA*R)/(c-c*BETA+c*BETA*n_m);















//GIVEN BY THE AUTHOR 

//h_squiggle = h;





//(CHECK IT ONCE AGAIN)

//epsilon_do_bar_squiggle = erfinv(-(erf(epsilon_do_bar/(SIGMA*2^0.5))+1.992))*SIGMA*2^0.5;
//f_do_bar_squiggle = 1/(2*SIGMA^2*pi)^0.5*exp(-epsilon_do_bar_squiggle^2/(2*SIGMA^2));
//F_do_bar_squiggle = (1 + erf(epsilon_do_bar_squiggle/(SIGMA*2^0.5)))/2;

//p_h = 1.4*y*12/h;
//n_do_squiggle = ((0.001-F_do_bar)*n_do*NU)/((1-NU)*(F_do_bar_squiggle-0.001));






end;

//calc. and check steady state
steady(maxit=100);

//check eigenvalues
check;

//standard deviation of shocks
shocks;
var epsilon_j = SIGMA^2;
var z_h = SIGMA^2;
end;

//stochastic simulation
stoch_simul(order=1,irf=40);
