% format long;
% clear;
% syms i m l b;
% 
% bet = .99; gam = 1.03; sig = .25; thet = .75; phit0 = 1.03; phit1 = 1.0;
% tau = (phit0/phit1) - 1; i_b = .2; g = .1; a = .8;
% 
% vars = [m b i l];
% 
% 
% q_l = bet*phit1*(m+(tau*m/sig)+l);
% q_b = bet*phit1*(m+(tau*m/sig)+b);
% l_q = a*inv(q_l+g)-1;
% S_B = (a*log(q_l+g) - q_l) - (a*log(q_b+g) - q_b) + phit1*bet*((1+i_b)*b - (1+i)*l);
% S_L = bet*phit1*i*(sig/(1-sig))*l;
% dl_di = l/(l_q-(1+i));
% di_dl = inv(dl_di);
% dpsi_di = bet*phit1*l_q*l*inv(l_q - (1+i));
% dpsi_dl = bet*phit1*l_q;
% 
% dSL_di = bet*phit1*(sig/(1-sig))*l;
% dSL_dL = bet*phit1*(sig/(1-sig))*i;
% dSB_di = - phit1*bet*l;
% dSB_dL = dpsi_dl - bet*phit1*(1+i);
% 
% nash_barg_l = -thet*S_L*dSB_dL == (1-thet)*S_B*dSL_dL;
% nash_barg_i = -thet*S_L*dSB_di == (1-thet)*S_B*dSL_di;
% % dm_priv = 1+i + di_dl*l == l_q;
% money_holdings = (gam-bet)/bet == sig*[l_q - 1];
% dm_pub = 1+i_b == a/(g+q_b);
% 
% kalai_barg = bet*phit1*(sig/(1-sig))*(l+(1+i)*dl_di) == 0;
% kalai_share = (1-thet)*S_L == thet*S_B;
% 
% 
% assume(i, 'positive'); assume(m, 'positive'); assume(l, 'positive'); 
% assume(b, 'positive'); assume(i <= i_b);
% 
% eqs_nash = [nash_barg_i, nash_barg_l, money_holdings, dm_pub];
% eqs_kalai = [kalai_barg, kalai_share, money_holdings, dm_pub];
% 
% 
% [i_star, l_star] = solve([nash_barg_l, nash_barg_i], [i, l],'PrincipalValue',true);
% 
% 
% [m_star, b_star, i_star, l_star] = solve(eqs_nash, vars, 'PrincipalValue', true)
% %[m_star, b_star, i_star, l_star] = solve(eqs_kalai, vars)
% q_star = bet*phit1*(m_star+(m_star*tau/sig) + l_star)

% 
% %% u(q) = x^alpha
% format long;
% clear;
% syms i m l b x;
% 
% bet = .99; gam = 1.03; sig = .25; thet = .5; phit0 = 1.03; phit1 = 1.0;
% tau = (phit0/phit1) - 1; i_b = .1; g = 5; a = .5;
% 
% vars = [m b i l];
% 
% q_fb = (1/(a*g))^(1/(a-1))
% d_star = q_fb/(phit1*bet)
% 
% q_l = bet*phit1*(m+(tau*m/sig)+l);
% q_b = bet*phit1*(m+(tau*m/sig)+b);
% u_q_l = g*q_l^a;
% u_q_b = g*q_b^a;
% l_q = a*((q_l)^(a-1))-1;
% l_q_b = a*((q_b)^(a-1))-1;
% 
% S_B = ((u_q_l - q_l) - phit1*bet*(1+i)*l)  - ((u_q_b - q_b) - phit1*bet*(1+i_b)*b);
% S_B2 = ((u_q_l - q_l) - phit1*bet*(l+x))  - ((u_q_b - q_b) - phit1*bet*(1+i_b)*b);
% S_L = bet*phit1*i*(sig/(1-sig))*l;
% S_L2 = bet*phit1*(sig/(1-sig))*x;
% dpsi_dl = bet*phit1*l_q;
% 
% dSL_di = bet*phit1*(sig/(1-sig))*l;
% dSL_dL = bet*phit1*(sig/(1-sig))*i;
% 
% dSL2_dL = 0;
% dSL2_dx = bet*phit1*(sig/(1-sig));
% 
% dSB_di = - phit1*bet*l;
% dSB_dL = dpsi_dl - bet*phit1*(1+i);
% 
% dSB2_dx = -bet*phit1;
% dSB2_dL = dpsi_dl - bet*phit1;
% 
% money_holdings = (gam-bet)/bet == sig*[l_q - 1];
% dm_pub = 1+i_b == l_q_b;
% 
% nash_barg_l = -thet*S_L*dSB_dL == (1-thet)*S_B*dSL_dL;
% nash_barg_i = -thet*S_L*dSB_di == (1-thet)*S_B*dSL_di;
% nash_barg_x = (1-thet)*(sig/(1-sig))*x+ thet*x == (thet/(bet*phit1))*(u_q_l - q_l - (u_q_b - q_b) + bet*phit1*((1+i_b)*b - l));
% nash_barg_l2 = -thet*S_L2*dSB2_dL == (1-thet)*S_B2*dSL2_dL;
% 
% kalai_barg = dpsi_dl - bet*phit1*(1+i) == 0;
% kalai_share = (1-thet)*S_L == thet*S_B;
% 
% % 
% % assume(i, 'positive'); assume(m, 'positive'); assume(l, 'positive'); 
% % assume(b, 'positive'); assume(i <= i_b);
% 
% eqs_nash = [nash_barg_i, nash_barg_l, money_holdings, dm_pub];
% eqs_nash_2 = [nash_barg_x, nash_barg_l2, money_holdings, dm_pub];
% eqs_kalai = [kalai_barg, kalai_share, money_holdings, dm_pub];
% 
% %[m_star, b_star, i_star, l_star] = vpasolve(eqs_nash, vars, [2,0.2,.05,0.2]) %vpa nash
% [m_star, b_star, i_star, l_star] = solve(eqs_nash, [m b i l]) %solve nash
% [m_star, b_star, x_star, l_star] = solve(eqs_nash_2, [m b x l]) %solve nash
% i_star = x_star/l_star;
% [m_star, b_star, i_star, l_star] = solve(eqs_kalai, vars) %vpa kalai
% 
% l_q_i = a*((bet*phit1*(m_star+(tau*m_star/sig)+l_star))^(a-1))-1;
% l_q_b = a*((bet*phit1*(m_star+(tau*m_star/sig)+b_star))^(a-1))-1;
% q_l = bet*phit1*(m_star+(m_star*tau/sig) + l_star);
% q_b = bet*phit1*(m_star+(m_star*tau/sig) + b_star);
% vpa([q_fb,0;
%     q_l, l_q_i;
%     q_b, l_q_b])

%% Corner Solution
format long;
clear;
syms i m l b x;

bet = .99; gam = 1.03; sig = .25; thet = .75; phit0 = 1.03; phit1 = 1.0;
tau = (phit0/phit1) - 1; i_b = .2; g = 2; a = .5;

vars = [m b i l];

q_fb = (1/(a*g))^(1/(a-1))
d_star = q_fb/(phit1*bet)

q_l = bet*phit1*(m+(tau*m/sig)+l);
q_b = bet*phit1*(m+(tau*m/sig)+b);
u_q_l = g*q_l^a;
u_q_b = g*q_b^a;
l_q = a*g*((q_l)^(a-1)) - 1;
l_q_b = a*g*((q_b)^(a-1)) - 1;

S_B = ((u_q_l - q_l) - phit1*bet*(1+i)*l)  - ((u_q_b - q_b) - phit1*bet*(1+i_b)*b);
S_B2 = ((u_q_l - q_l) - phit1*bet*(l+x))  - ((u_q_b - q_b) - phit1*bet*(1+i_b)*b);
S_L = bet*phit1*i*(sig/(1-sig))*l;
S_L2 = bet*phit1*(sig/(1-sig))*x;
dpsi_dl = bet*phit1*l_q;

dSL_di = bet*phit1*(sig/(1-sig))*l;
dSL_dL = bet*phit1*(sig/(1-sig))*i;

dSL2_dL = 0;
dSL2_dx = bet*phit1*(sig/(1-sig));

dSB_di = - phit1*bet*l;
dSB_dL = dpsi_dl - bet*phit1*(1+i);

dSB2_dx = -bet*phit1;
dSB2_dL = dpsi_dl - bet*phit1;

money_holdings = (gam-bet)/bet == sig*[l_q - 1];
dm_pub = 1+i_b == l_q_b;

nash_barg_l = -thet*S_L*dSB_dL == (1-thet)*S_B*dSL_dL;
nash_barg_i = -thet*S_L*dSB_di == (1-thet)*S_B*dSL_di;

nash_barg_x = x == (thet/(bet*phit1*((1-thet)*(sig/(1-sig))+ thet)))*(u_q_l - q_l - (u_q_b - q_b) + bet*phit1*((1+i_b)*b - l));
nash_barg_l2 = l == ((2/(a*g))^(1/(a-1)))/(bet*phit1) - (1+tau/sig)*m;

kalai_barg = dpsi_dl - bet*phit1*(1+i) == 0;
kalai_share = (1-thet)*S_L == thet*S_B;

% 
% assume(i, 'positive'); assume(m, 'positive'); assume(l, 'positive'); 
% assume(b, 'positive'); assume(i <= i_b);

eqs_nash = [nash_barg_i, nash_barg_l, money_holdings, dm_pub];
eqs_nash_2 = [nash_barg_l2, nash_barg_x, money_holdings, dm_pub];
eqs_kalai = [kalai_barg, kalai_share, money_holdings, dm_pub];

%[m_star, b_star, i_star, l_star] = vpasolve(eqs_nash, vars, [2,0.2,.05,0.2]) %vpa nash
[m_star, b_star, i_star, l_star] = solve(eqs_nash, [m b i l]) %solve nash
[m_star, b_star, x_star, l_star] = solve(eqs_nash_2, [m b x l]) %solve nash
i_star = x_star/l_star;
[m_star, b_star, i_star, l_star] = solve(eqs_kalai, vars) %vpa kalai

l_q_i = a*((bet*phit1*(m_star+(tau*m_star/sig)+l_star))^(a-1))-1;
l_q_b = a*((bet*phit1*(m_star+(tau*m_star/sig)+b_star))^(a-1))-1;
q_l = bet*phit1*(m_star+(m_star*tau/sig) + l_star);
q_b = bet*phit1*(m_star+(m_star*tau/sig) + b_star);
vpa([q_fb,0;
    q_l, l_q_i;
    q_b, l_q_b])