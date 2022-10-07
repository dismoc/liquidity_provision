%% DW Rate
syms x;
format long;
clear;

bin = 100;

bet = .99; sig = .21; thet = .52; phit0 = 1.02; phit1 = 1;
tau = (phit0/phit1) - 1; i_b = [0:.1/bin:.1]; g =1; a = .6; gam = 1 + tau;

q_star = (a*g)^(1-a); 
q_opt = ((2*a*sig*bet)/(gam-bet*(1-2*sig)))^(1-a);
%q_opt_ns2 = (((gam-bet)/(bet*sig) + 1)/(g*a))^(1/(a-1))
m_opt = q_opt/(2*bet*phit1);

k=1;

for k = 1:length(i_b)
    m(k) = q_opt/(2*bet*phit1);
    %m_ns2(k) = q_opt_ns2*bet*phit1;
    
    l(k) = m(k);
    d_b(k) = (((a*g)/(1+i_b(k)))^(1-a))/(phit1*bet);
    b_res(k) = max(d_b(k) - m(k),0);
    
    q_l(k) = bet*phit1*(m(k)+l(k));
    q_b(k) = bet*phit1*(m(k)+b_res(k));
    
    u = @(q) g*q^a;
    psi_l(k) = u(q_l(k)) - q_l(k);
    s_b(k) = u(q_b(k)) - q_b(k) - bet*phit1*(i_b(k))*b_res(k);
    x(k) = max((thet/bet*phit1)*(psi_l(k) - s_b(k)),0);
    i(k) = x(k)/l(k);
    
    %Welfare
    if 2*m(k) < d_b
        m(k) = 0; l(k) = 0; x(k) = 0; i(k) = 0;
        b(k) = d_b(k) - m(k) - l(k);
        q_l_b(k) = bet*phit1*(m(k) + l(k) + b(k));
        W_b(k) = sig*(u(q_b(k)) - q_b(k));
    else
        b(k) = 0;
        q_l_b(k) = 0;
        W_b(k) = sig*(u(q_b(k)) - q_b(k));
    end
    W_l(k) = sig*(u(q_l(k)) - q_l(k));
    W(k) = max(W_l(k),W_b(k));
    %W_ns2(k) = (u(m_ns2(k))-q_opt_ns2)*sig ;
end



close;
subplot(2,3,1); plot(i_b,W,i_b,sig*x,'--'); ylabel('Welfare'); xlabel('DW Rate');
subplot(2,3,3); plot(i_b,m); ylabel('Money holdings'); xlabel('DW Rate');
subplot(2,3,6); plot(i_b,i./i_b); ylabel('i^l as % of i^b'); xlabel('DW Rate');
subplot(2,3,2); plot(i_b,x); ylabel('Transfer to Lender - x'); xlabel('DW Rate');
subplot(2,3,5); plot(i_b,i); ylabel('Private Interest - i^l'); xlabel('DW Rate');
subplot(2,3,4); plot(i_b,b); ylabel('Public Borrowing - b'); xlabel('DW Rate');

print -djpeg epsFig1case
%% Welfare Diff alpha
clear;

bin = 100;

bet = .99; sig = .21; thet = .52; phit0 = 1.02; phit1 = 1;
tau = (phit0/phit1) - 1; i_b = .0075; g =1; a = [0.01:.01:.99] ; gam = 1 + tau;

for k = 1:length(a)
    
    q_opt = ((1 + (gam-bet*(1-2*sig))/(bet*sig*2))/(g*a(k)))^(1/(a(k)-1));
    q_opt_ns2 = (((gam-bet)/(bet*sig) + 1)/(g*a(k)))^(1/(a(k)-1))
    m_opt = q_opt/(2*bet*phit1);

    % money holdings
    m(k) = q_opt/(2*bet*phit1);
    m_ns2(k) = q_opt_ns2*bet*phit1;
    l(k) = m(k);
    
    d_b(k) = (((2 + i_b)/(a(k)*g))^(1/(a(k)-1)))/(phit1*bet);
    b_res(k) = max(d_b(k) - m(k),0);
    
    q_l(k) = bet*phit1*(m(k)+l(k));
    q_b(k) = bet*phit1*(m(k)+b_res(k));
    
    u = @(q) g*q^a(k);
    
    psi_l(k) = u(q_l(k)) - q_l(k);
    s_b(k) = u(q_b(k)) - q_b(k) - bet*phit1*(1+i_b)*b_res(k);
    x(k) = max((thet/bet*phit1)*(psi_l(k) - bet*phit1*l(k) - s_b(k)),0);
    i(k) = x(k)/l(k);
    
    
    %Welfare
    if 2*m(k) < d_b
        b(k) = d_b(k) - m(k) - l(k);
        q_l_b(k) = bet*phit1*(m(k) + l(k) + b(k));
        W_b(k) = sig*(u(q_l_b(k)) - q_l_b(k));
    else
        b(k) = 0;
        q_l_b(k) = 0;
        W_b(k) = sig*(u(q_b(k)) - q_b(k));
    end
    W_l(k) = sig*(u(q_l(k)) - q_l(k));
    W(k) = max(W_l(k),W_b(k));
    W_ns2(k) = (u(m_ns2(k))-q_opt_ns2)*sig;
end

close;
subplot(2,1,1); plot(a,W - W_ns2,a,zeros(1,length(a)),':'); xlabel('Degree of Homogeneity'); ylabel('Diff in Welfare(Open - Close)');
subplot(2,1,2); plot(a,W,a,W_ns2,'--'); 
xlabel('Degree of Homogeneity'); ylabel('Welfare');
legend('Open','Close');

print -djpeg epsFig2
%% Welfare Diff sigma
clear;

bin = 100;

bet = .99; sig = [0:.005:.5]; thet = .52; phit0 = 1.02; phit1 = 1;
tau = (phit0/phit1) - 1; i_b = .02; g =1; a = .33 ; gam = 1 + tau;

for k = 1:length(sig)
    
    q_opt = ((1 + (gam-bet*(1-2*sig(k)))/(bet*sig(k)*2))/(g*a))^(1/(a-1));
    q_opt_ns2 = (((gam-bet)/(bet*sig(k)) + 1)/(g*a))^(1/(a-1))
    m_opt = q_opt/(2*bet*phit1);

    % money holdings
    m(k) = q_opt/(2*bet*phit1);
    m_ns2(k) = q_opt_ns2*bet*phit1;
    l(k) = m(k);
    
    d_b(k) = (((2 + i_b)/(a*g))^(1/(a-1)))/(phit1*bet);
    b_res(k) = max(d_b(k) - m(k),0);
    
    q_l(k) = bet*phit1*(m(k)+l(k));
    q_b(k) = bet*phit1*(m(k)+b_res(k));
    
    u = @(q) g*q^a;
    
    psi_l(k) = u(q_l(k)) - q_l(k);
    s_b(k) = u(q_b(k)) - q_b(k) - bet*phit1*(1+i_b)*b_res(k);
    x(k) = max((thet/bet*phit1)*(psi_l(k) - bet*phit1*l(k) - s_b(k)),0);
    i(k) = x(k)/l(k);
    
    
    %Welfare
    if 2*m(k) < d_b
        b(k) = d_b(k) - m(k) - l(k);
        q_l_b(k) = bet*phit1*(m(k) + l(k) + b(k));
        W_b(k) = sig(k)*(u(q_l_b(k)) - q_l_b(k));
    else
        b(k) = 0;
        q_l_b(k) = 0;
        W_b(k) = sig(k)*(u(q_b(k)) - q_b(k));
    end
    W_l(k) = sig(k)*(u(q_l(k)) - q_l(k));
    W(k) = max(W_l(k),W_b(k));
    W_ns2(k) = (u(m_ns2(k))-q_opt_ns2)*sig(k);
end


close;
subplot(2,1,1); plot(sig,W - W_ns2,sig,zeros(1,length(sig)),':'); xlabel('Meeting Prob'); ylabel('Diff in Welfare(Open - Close)');
subplot(2,1,2); plot(sig,W,sig,W_ns2,'--'); 
xlabel('Meeting Prob'); ylabel('Welfare');
legend('Open','Close');


print -djpeg epsFig5
