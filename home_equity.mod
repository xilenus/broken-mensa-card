//endogenous variables
var
c h d 
//c_squiggle h_squiggle d_squiggle

epsilon_bar epsilon_do_bar
//epsilon_bar_squiggle epsilon_do_bar_squiggle

n_m n_u n_s n_do
//n_m_squiggle n_u_squiggle n_s_squiggle n_do_squiggle

n_u_hat n_s_hat m_hat g_hat g_f_hat v_hat

//V
p_h R
lambda_cc

f_do_bar F_do_bar f_bar F_bar
//f_do_bar_squiggle F_do_bar_squiggle f_bar_squiggle F_bar_squiggle
f_j F_j

//lambda_t_squiggle
y z_a

G //G_squiggle
;

//exogenous variables
varexo
epsilon_j z_h
;

//parameters
parameters
BETA GAMMA ALPHA_imp ALPHA_pat PSI SIGMA ZETA RHO_u ETA MU THETA OMEGA XI KAPPA NU CHI PHI_h PHI_a pi expon





z_a_ST
p_h_ST
h_ST

epsilon_do_bar_ST
F_do_bar_ST
epsilon_bar_ST
F_bar_ST
epsilon_j_ST
F_j_ST
n_s_ST
n_u_ST
g_hat_ST
G_ST
y_ST
c_ST
d_ST
lambda_cc_ST
R_ST
n_m_ST


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
//THETA = 0.181;
THETA = 0;
OMEGA = 0.3333;
//XI = 0.98;
XI = 1;
KAPPA = -4.381;
//NU = 0.2;
NU = 1;
CHI = 0.8;
PHI_h = 0.983;
PHI_a = 0.983;




















//model

model;

//optimization problem for impatient households

//f_do_bar = 1/(2*SIGMA^2*pi)^0.5*expon^(-epsilon_do_bar^2/(2*SIGMA^2));
//f_do_bar_ST*f_do_bar = -epsilon_do_bar_ST^2*epsilon_do_bar*exp(-(epsilon_do_bar_ST^2)/(2*SIGMA^2))/(2*pi*SIGMA^6)^0.5;
epsilon_do_bar_ST = -SIGMA*(f_do_bar/epsilon_do_bar)^0.5;

//F_do_bar = (1 + erf(epsilon_do_bar/(SIGMA*2^0.5)))/2;
//F_do_bar_ST*F_do_bar = epsilon_do_bar_ST*epsilon_do_bar*exp(-(epsilon_do_bar_ST^2)/(SIGMA^2*2))/(SIGMA*2^0.5*pi^0.5);
F_do_bar_ST = (epsilon_do_bar_ST*epsilon_do_bar*exp(-(epsilon_do_bar_ST^2)/(SIGMA^2*2))/(SIGMA*2^0.5*pi^0.5))/F_do_bar;

//f_bar = 1/(2*SIGMA^2*pi)^0.5*expon^(-epsilon_bar^2/(2*SIGMA^2));
//f_bar_ST*f_bar = -epsilon_bar_ST^2*epsilon_bar*exp(-(epsilon_bar_ST^2)/(2*SIGMA^2))/(2*pi*SIGMA^6)^0.5;
epsilon_bar_ST = -SIGMA*(f_bar/epsilon_bar)^0.5;

//F_bar = (1 + erf(epsilon_bar/(SIGMA*2^0.5)))/2;
//F_bar_ST*F_bar = epsilon_bar_ST*epsilon_bar*exp(-(epsilon_bar_ST^2)/(SIGMA^2*2))/(SIGMA*2^0.5*pi^0.5);
F_bar_ST = (epsilon_bar_ST*epsilon_bar*exp(-(epsilon_bar_ST^2)/(SIGMA^2*2))/(SIGMA*2^0.5*pi^0.5))/F_bar;

//f_j = 1/(2*SIGMA^2*pi)^0.5*expon^(-epsilon_j^2/(2*SIGMA^2));
//f_j_ST*f_j = -epsilon_j_ST^2*epsilon_j*exp(-(epsilon_j_ST^2)/(2*SIGMA^2))/(2*pi*SIGMA^6)^0.5;
//epsilon_j_ST = -SIGMA*(f_j/epsilon_j)^0.5;
epsilon_j_ST = 0;

//F_j = (1 + erf(epsilon_j/(SIGMA*2^0.5)))/2;
//F_j_ST*F_j = epsilon_j_ST*epsilon_j*exp(-(epsilon_j_ST^2)/(SIGMA^2*2))/(SIGMA*2^0.5*pi^0.5);
//F_j_ST = (epsilon_j_ST*epsilon_j*exp(-(epsilon_j_ST^2)/(SIGMA^2*2))/(SIGMA*2^0.5*pi^0.5))/F_j;
F_j_ST = 0;



//n_s = n_u + RHO_u*(1-n_u);
//n_s_ST*n_s = n_u_ST*n_u*(1-RHO_u);
n_s_ST = -n_u*RHO_u/(n_s-n_u);

//n_u = n_s(-1)*(1-g_hat(-1)+OMEGA*g_hat(-1)*(1-F_do_bar));
//n_u*n_u_ST = n_s_ST*n_s(-1)*(1-g_hat_ST+OMEGA*g_hat_ST*(1-F_do_bar_ST)) + F_do_bar_ST*F_do_bar*(-OMEGA*g_hat_ST*n_s_ST) + g_hat_ST*g_hat(-1)*(-n_s_ST + n_s_ST*OMEGA*(1-F_do_bar_ST));
//n_u_ST = (n_s_ST^2*g_hat(-1)-F_do_bar_ST*n_s_ST*g_hat_ST*OMEGA*F_do_bar)/(n_u-n_s(-1)-n_u_ST*g_hat(-1));
n_u_ST = -(n_s_ST^2*g_hat(-1)+OMEGA*n_s_ST*g_hat_ST*F_do_bar_ST*F_do_bar)/(n_u-n_s(-1)-n_s_ST*g_hat(-1));

//n_m = OMEGA*g_hat(-1)*n_s(-1)*F_do_bar + (1-OMEGA*g_hat(-1)*n_s(-1))*F_bar;
//n_m_ST*n_m = F_do_bar_ST*F_do_bar*OMEGA*g_hat_ST*n_s_ST + F_bar_ST*F_bar*(1-OMEGA*g_hat_ST*n_s_ST) + g_hat_ST*g_hat(-1)*(OMEGA*n_s_ST*F_do_bar_ST-OMEGA*F_bar_ST*n_s_ST) + n_s_ST*n_s(-1)*(OMEGA*g_hat_ST*F_do_bar_ST-F_bar_ST*OMEGA*g_hat_ST);
g_hat_ST = (n_m_ST*n_m-n_m_ST*F_bar-(F_do_bar_ST-F_bar_ST)*(g_hat(+1)+1)*(n_m_ST-1)*F_bar_ST/F_do_bar_ST)/(OMEGA*((F_bar_ST*(F_do_bar_ST-F_bar_ST)*(g_hat(-1)+1))/F_do_bar_ST-n_s*F_do_bar_ST));

//OTHER EQUATIONS


//G = BETA*(-OMEGA*g_hat*(epsilon_do_bar(+1)^2*(1+erf(epsilon_do_bar(+1)/(SIGMA*2^0.5)))/4 - epsilon_bar(+1)^2*(1+erf(epsilon_bar(+1)/(SIGMA*2^0.5)))/4)+(epsilon_do_bar(+1)-epsilon_bar(+1))*(1-g_hat+OMEGA*g_hat*(1-F_do_bar(+1)))+epsilon_bar(+1)*OMEGA*g_hat*(F_do_bar(+1)-F_bar(+1)));
G_ST*G/BETA = g_hat_ST*g_hat*(-OMEGA*(epsilon_do_bar_ST^2*(1+erf(epsilon_do_bar_ST/(SIGMA*2^0.5)))/4-epsilon_bar_ST^2*(epsilon_bar_ST/(SIGMA*2^0.5))/4)+(OMEGA*(1-F_do_bar_ST)-1)*(epsilon_do_bar_ST-epsilon_bar_ST)) + F_bar_ST*F_bar(+1)*(-OMEGA*g_hat_ST*epsilon_bar_ST) + F_do_bar_ST*F_do_bar(+1)*(-OMEGA*g_hat_ST+epsilon_bar_ST*OMEGA*g_hat_ST) + epsilon_do_bar_ST*epsilon_do_bar(+1)*(-OMEGA*g_hat_ST*(1/2*epsilon_do_bar_ST*(1+erf(epsilon_do_bar_ST/(SIGMA*2^0.5)))+epsilon_do_bar_ST^2*exp(-epsilon_do_bar_ST^2/(2*SIGMA^2))/(2*2^0.5*pi^0.5*SIGMA))+1-g_hat_ST+OMEGA*g_hat_ST*(1-F_do_bar_ST)) + epsilon_bar_ST*epsilon_bar(+1)*(-OMEGA*g_hat_ST*(1/2*epsilon_bar_ST*(erf(epsilon_bar_ST/(SIGMA*2^0.5))+1)+epsilon_bar_ST^2*exp(-epsilon_bar_ST^2/(2*SIGMA^2))/(2*2^0.5*pi^0.5*SIGMA))-(1-g_hat_ST+OMEGA*g_hat_ST*(1-F_do_bar_ST)));




y_ST = XI*z_a_ST;
c_ST = y_ST/(epsilon_do_bar_ST-epsilon_bar_ST+KAPPA+(1-RHO_u)*G_ST);
d_ST = n_m_ST*p_h_ST*h_ST*CHI/(2-n_m_ST);
lambda_cc_ST = (ZETA/c_ST+epsilon_bar_ST-PSI)/(CHI*p_h_ST*h_ST-d_ST);
R_ST = (1/c_ST-lambda_cc_ST+BETA*lambda_cc_ST*(1-n_m_ST))*c_ST/BETA;

//d = n_m*CHI*p_h*h + (1-n_m)*d(-1);                                          p_h
//d_ST*d = n_m_ST*n_m*CHI*p_h_ST*h_ST + p_h_ST*p_h*n_m_ST*CHI*h_ST + h_ST*h*n_m_ST*CHI*p_h_ST - n_m_ST*n_m*d_ST + d_ST*d(-1)*(1-n_m_ST);
n_m_ST = (d - d(-1))/(p_h+h-1);






//y = XI*z_a;
//y_ST*y = z_a_ST*z_a*XI;
y = z_a;


//c + p_h*(h-h(-1)) + ZETA*n_m + R(-1)*d(-1) = (1-n_u)*z_a*XI + d;          c
c_ST*c + h_ST*h*p_h_ST - h_ST*h*p_h_ST - h_ST*h(-1)*p_h_ST + n_m*n_m_ST*ZETA + R_ST*R(-1)*d_ST + d_ST*d(-1)*R_ST = d_ST*d + z_a_ST*z_a*(1-n_u_ST)*XI - n_u_ST*n_u*z_a_ST*XI;



//ALPHA_imp*z_h/h - p_h/c + BETA*p_h(+1)/c(+1) + lambda_cc*n_m*CHI*p_h = 0;   lambda_cc
//z_h_ST*z_h*ALPHA_imp/h_ST - h_ST*h*ALPHA_imp*z_h_ST/h_ST^2 + p_h_ST*p_h*(lambda_cc_ST*n_m_ST*CHI - 1/c_ST) + c_ST*c*p_h_ST/c_ST^2 + p_h_ST*p_h(+1)*BETA/c_ST - c_ST*c(+1)*BETA*p_h_ST/c_ST^2 + lambda_cc_ST*lambda_cc*n_m_ST*CHI*p_h_ST + n_m_ST*n_m*lambda_cc_ST*CHI*p_h_ST = 0;
p_h_ST*p_h*(lambda_cc_ST*n_m_ST*CHI - 1/c_ST) + c_ST*c*p_h_ST/c_ST^2 + p_h_ST*p_h(+1)*BETA/c_ST - c_ST*c(+1)*BETA*p_h_ST/c_ST^2 + lambda_cc_ST*lambda_cc*n_m_ST*CHI*p_h_ST + n_m_ST*n_m*lambda_cc_ST*CHI*p_h_ST = 0;


//1/c - lambda_cc + BETA*(lambda_cc(+1)*(1-n_m(+1))-R/c(+1)) = 0;             R
-c/c_ST - lambda_cc_ST*lambda_cc + lambda_cc_ST*lambda_cc(+1)*(1-n_m_ST)*BETA - n_m_ST*n_m(+1)*lambda_cc_ST*BETA - R_ST*R/c_ST + c_ST*c(+1)*R_ST/c_ST^2 = 0;

//PSI = ZETA/c + epsilon_bar - lambda_cc*(CHI*p_h*h-d(-1));                   epsilon_bar
c_ST*c*(-ZETA/c_ST^2) + epsilon_bar_ST*epsilon_bar + lambda_cc_ST*lambda_cc*(d_ST-CHI*p_h_ST*h_ST) - p_h_ST*p_h*lambda_cc_ST*CHI*h_ST - h_ST*h*lambda_cc_ST*CHI*p_h_ST + d_ST*d(-1)*lambda_cc_ST = 0;

//epsilon_do_bar - epsilon_bar = y/c - KAPPA + (1-RHO_u)*G;                   epsilon_do_bar
epsilon_do_bar_ST*epsilon_do_bar - epsilon_bar_ST*epsilon_bar = y_ST*y/c_ST - c*y_ST/c_ST + G_ST*G*(1-RHO_u);



//n_do = OMEGA*g_hat(-1)*n_s;
//n_do_ST*n_do = g_hat_ST*g_hat(-1)*OMEGA*n_s_ST + n_s_ST*n_s*OMEGA*g_hat_ST;
n_do = g_hat(-1)+n_s;

//LABOR MARKET

//n_u_hat = NU*n_u;
//n_u_hat_ST*n_u_hat = n_u_ST*n_u*NU;
n_u_hat = n_u;

//n_s_hat = NU*n_s;
//n_s_hat_ST*n_s_hat = NU*n_s_ST*n_s;
n_s_hat = n_s;

//m_hat = MU*n_s_hat^ETA*v_hat^(1-ETA);
//m_hat_ST*m_hat = n_s_hat*ETA*MU*v_hat_ST^(1-ETA)*n_s_hat_ST^ETA;
m_hat = n_s_hat*ETA;

//g_hat = m_hat/n_s_hat;
//g_hat_ST*g_hat = m_hat_ST*m_hat/n_s_hat_ST + n_s_hat_ST*n_s_hat*(-m_hat_ST/n_s_hat_ST^2);
//g_hat = m_hat - n_s_hat;
//g_f_hat = m_hat/v_hat;                  v_hat
//g_f_hat_ST*g_f_hat = m_hat_ST*m_hat/v_hat_ST - v_hat*m_hat_ST/v_hat_ST;
//g_f_hat = m_hat - v_hat;
g_f_hat = g_hat + n_s_hat - v_hat;




//FIRMS

//0 = (1 - XI)*z_a;                                               z_a
//0 = z_a_ST*z_a*(1-XI);
//n_s_hat = n_s;

//THETA = g_f_hat*(1-OMEGA+OMEGA*NU*n_s/n_s_hat*F_do_bar(+1))*0;  g_f_hat
//0 = g_f_hat_ST*g_f_hat*(1-OMEGA+OMEGA*NU*n_s_ST/n_s_hat_ST*F_do_bar_ST)*0 + (n_s_ST*n_s*OMEGA*NU/n_s_hat_ST*F_do_bar_ST*0 + n_s_hat*(-OMEGA*NU*n_s/n_s_hat_ST*f_do_bar_ST)*0+F_do_bar_ST*F_do_bar(+1)*OMEGA*NU*n_s/n_s_hat_ST)*g_f_hat_ST;


//NU*h = 1;                                                       h
//h_ST*h*NU = 0;
//h = 0;


//NU*d = 1000;                                                    d
//d_ST*d*NU = 0;
//d = 0;



end;






//epsilon_do_bar_ST = -SIGMA*(f_do_bar/epsilon_do_bar)^0.5;
epsilon_do_bar_ST = -SIGMA;

//F_do_bar_ST = (epsilon_do_bar_ST*epsilon_do_bar*exp(-(epsilon_do_bar_ST^2)/(SIGMA^2*2))/(SIGMA*2^0.5*pi^0.5))/F_do_bar;
F_do_bar_ST = 0;

//epsilon_bar_ST = -SIGMA*(f_bar/epsilon_bar)^0.5;
epsilon_bar_ST = -SIGMA;

//F_bar_ST = (epsilon_bar_ST*epsilon_bar*exp(-(epsilon_bar_ST^2)/(SIGMA^2*2))/(SIGMA*2^0.5*pi^0.5))/F_bar;
F_bar_ST = 0;

epsilon_j_st = 0;
F_j_ST = 0;

//g_hat_ST = (n_m_ST*n_m-n_m_ST*F_bar-(F_do_bar_ST-F_bar_ST)*(g_hat(+1)+1)*(n_m_ST-1)*F_bar_ST/F_do_bar_ST)/(OMEGA*((F_bar_ST*(F_do_bar_ST-F_bar_ST)*(g_hat(-1)+1))/F_do_bar_ST-n_s*F_do_bar_ST));
g_hat_ST = (-(F_do_bar_ST-F_bar_ST)*(n_m_ST-1)*F_bar_ST/F_do_bar_ST)/(OMEGA*((F_bar_ST*(F_do_bar_ST-F_bar_ST))/F_do_bar_ST));

//n_s_ST = -n_u*RHO_u/(n_s-n_u);
n_s_ST  = -RHO_u;

//n_u_ST = -(n_s_ST^2*g_hat(-1)+OMEGA*n_s_ST*g_hat_ST*F_do_bar_ST*F_do_bar)/(n_u-n_s(-1)-n_s_ST*g_hat(-1));
n_u_ST = n_s_ST+OMEGA*g_hat_ST*F_do_bar_ST;



G_ST = BETA*(((g_hat_ST*(-OMEGA*(epsilon_do_bar_ST^2*(1+erf(epsilon_do_bar_ST/(SIGMA*2^0.5)))/4-epsilon_bar_ST^2*(epsilon_bar_ST/(SIGMA*2^0.5))/4)+(OMEGA*(1-F_do_bar_ST)-1)*(epsilon_do_bar_ST-epsilon_bar_ST)) + F_bar_ST*(-OMEGA*g_hat_ST*epsilon_bar_ST)) + epsilon_do_bar_ST*(-OMEGA*g_hat_ST*(1/2*epsilon_do_bar_ST*(1+erf(epsilon_do_bar_ST/(SIGMA*2^0.5)))+epsilon_do_bar_ST^2*exp(-epsilon_do_bar_ST^2/(2*SIGMA^2))/(2*2^0.5*pi^0.5*SIGMA))+1-g_hat_ST+OMEGA*g_hat_ST*(1-F_do_bar_ST)) + epsilon_bar_ST*(-OMEGA*g_hat_ST*(1/2*epsilon_bar_ST*(erf(epsilon_bar_ST/(SIGMA*2^0.5))+1)+epsilon_bar_ST^2*exp(-epsilon_bar_ST^2/(2*SIGMA^2))/(2*2^0.5*pi^0.5*SIGMA))-(1-g_hat_ST+OMEGA*g_hat_ST*(1-F_do_bar_ST))) + F_do_bar_ST*(-OMEGA*g_hat_ST+epsilon_bar_ST*OMEGA*g_hat_ST)));

y_ST = XI*z_a_ST;
c_ST = y_ST/(epsilon_do_bar_ST-epsilon_bar_ST+KAPPA+(1-RHO_u)*G_ST);
d_ST = n_m_ST*p_h_ST*h_ST*CHI/(2-n_m_ST);
lambda_cc_ST = (ZETA/c_ST+epsilon_bar_ST-PSI)/(CHI*p_h_ST*h_ST-d_ST);
R_ST = (1/c_ST-lambda_cc_ST+BETA*lambda_cc_ST*(1-n_m_ST))*c_ST/BETA;

n_m_ST = 0;






//assign steady state values
initval;




c = 0;
h = 0;
d = 0;
epsilon_bar = 0;
epsilon_do_bar = 0;
n_m = 0;
n_u = 0;
n_s = 0;
n_do = 0;
n_u_hat = 0;
n_s_hat = 0;
m_hat = 0;
g_hat = 0;
g_f_hat = 0;
v_hat = 0;
p_h = 0;
R = 0;
lambda_cc = 0;
f_do_bar = 0;
F_do_bar = 0;
f_bar = 0;
F_bar = 0;
f_j = 0;
F_j = 0;
y = 0;
z_a = 0;
G = 0;













end;

//calc. and check steady state
steady(maxit=10000);

//check eigenvalues
check;

//standard deviation of shocks
shocks;
var epsilon_j = SIGMA^2;
var z_h = SIGMA^2;
end;

//stochastic simulation
stoch_simul(order=1,irf=40);
