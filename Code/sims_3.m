%% DW Rate
syms x;
format long;
clear;

bin = 100;

bet = .98; sig = [.01:.01:.05]; thet = .8; phit0 = 1.02; phit1 = 1;
tau = (phit0/phit1) - 1; i_b = [0:.1/bin:.2]; g =2; a = .8; gam = 1 + tau;

q_opt = ((1 + gam/bet)/(g*a))^(1/(a-1));
q_opt_ns2 = (((gam-bet)/(bet*sig) + 1)/(g*a))^(1/(a-1))
m_opt = q_opt*sig/(bet*phit1);

for i=1:length(sig)
    y(i) = (gam-bet*(1-2*sig(i)))/(2*bet*sig(i));
    y2(i) = gam/bet - 1;
end
plot(sig,y,sig,y2)

for k = 1:length(i_b)
    m(k) = q_opt*sig/(bet*phit1);
    m_ns2(k) = q_opt_ns2*bet*phit1;
    l(k) = (1-sig)*m(k)/sig;
    d_b(k) = (((2 + i_b(k))/(a*g))^(1/(a-1)))/(phit1*bet);
    b_res(k) = max(d_b(k) - m(k),0);
    q_l(k) = bet*phit1*m(k)/sig;
    q_b(k) = bet*phit1*(m(k)+b_res(k));
    u = @(q) g*q^a;
    psi_m(k) = u(q_l(k)) - q_l(k);
    s_b(k) = u(q_b(k)) - q_b(k) - bet*phit1*(1+i_b(k))*b_res(k);
    x(k) = max((thet/bet*phit1)*(psi_m(k) - bet*phit1*l(k) - s_b(k)),0);
    i(k) = x(k)/l(k);
    if m(k)/sig < d_b
        b(k) = d_b(k) - m(k)/sig;
        q_l_b(k) = bet*phit1*(m(k)/sig + b(k));
        W_b(k) = u(q_l_b(k)) - q_l_b(k) - phit0*m(k);
    else
        b(k) = 0;
        q_l_b(k) = 0;
        W_b(k) = u(q_b(k)) - q_b(k) - phit0*m(k);
    end
    W_l(k) = u(q_l(k)) - q_l(k) - phit0*m(k);
    W(k) = max(W_l(k),W_b(k));
    W_ns2(k) = u(m_ns2(k))*sig + (1-sig)*(phit1-phit0)*m_ns2(k);
    
    if k > 1
        cp_elas(k-1) = ((i(k) - i(k-1))/((i(k) + i(k-1))/2))/((i_b(k) - i_b(k-1))/((i_b(k) + i_b(k-1))/2));
    end
end
        cp_elas(k) = ((i(k) - i(k-1))/((i(k) + i(k-1))/2))/((i_b(k) - i_b(k-1))/((i_b(k) + i_b(k-1))/2));

close; sgtitle('Effect of monetary policy');
subplot(2,3,1); plot(i_b,sig*W,i_b,x,'--'); ylabel('Welfare'); xlabel('DW Rate');
subplot(2,3,6); plot(i_b,i./i_b); ylabel('i^l as % of i^b'); xlabel('DW Rate');
subplot(2,3,3); plot(i_b,m); ylabel('Money Holdings - m '); xlabel('DW Rate');
subplot(2,3,2); plot(i_b,x); ylabel('Transfer to Lender - x'); xlabel('DW Rate');
subplot(2,3,5); plot(i_b,i); ylabel('Private Interest - i^l'); xlabel('DW Rate');
subplot(2,3,4); plot(i_b,b); ylabel('Public Borrowing - b'); xlabel('DW Rate');

print -deps epsFig1
%% Welfare Diff alpha
clear;

bin = 100;

bet = .98; sig = .1; thet = .8; phit0 = 1.02; phit1 = 1;
tau = (phit0/phit1) - 1; i_b = .02; g =1.4; a = [0.01:.01:.99] ; gam = 1 + tau;

for k = 1:length(a)
    
    q_opt = ((1 + gam/bet)/(g*a(k)))^(1/(a(k)-1));
    q_opt_ns2 = (((gam-bet)/(bet*sig) + 1)/(g*a(k)))^(1/(a(k)-1));
    m_opt = q_opt*sig/(bet*phit1);

    m(k) = q_opt*sig/(bet*phit1);
    m_ns2(k) = q_opt_ns2*bet*phit1;
    
    l(k) = (1-sig)*m(k)/sig;
    d_b(k) = (((2 + i_b)/(a(k)*g))^(1/(a(k)-1)))/(phit1*bet);
    b_res(k) = max(d_b(k) - m(k),0);
    q_l(k) = bet*phit1*m(k)/sig;
    q_b(k) = bet*phit1*(m(k)+b_res(k));
    u = @(q) g*q^a(k);
    psi_m(k) = u(q_l(k)) - q_l(k);
    s_b(k) = u(q_b(k)) - q_b(k) - bet*phit1*(1+i_b)*b_res(k);
    x(k) = max((thet/bet*phit1)*(psi_m(k) - bet*phit1*l(k) - s_b(k)),0);
    i(k) = x(k)/l(k);
    
    
    if m(k)/sig < d_b
        b(k) = d_b(k) - m(k)/sig;
        q_l_b(k) = bet*phit1*(m(k)/sig + b(k));
        W_b(k) = sig*(u(q_l_b(k)) - q_l_b(k)) - (gam-bet)*m(k);
    else
        b(k) = 0;
        q_l_b(k) = 0;
        W_b(k) = sig*(u(q_b(k)) - q_b(k)) - (gam-bet)*m(k);
    end
    W_l(k) = sig*(u(q_l(k)) - q_l(k)) - (gam-bet)*m(k);
    W(k) = max(W_l(k),W_b(k));
    
    q_ns2(k) = phit1*bet*m_ns2(k);
    W_ns2(k) = (u(q_ns2(k))-q_ns2(k))*sig - (gam-bet)*m_ns2(k);
    
    q_fb(k) = inv(a(k)*g)^(inv(a(k)-1));
    m_fb(k) = q_fb(k)/(bet*phit1);
    W_fb(k) = sig*(u(q_fb(k)) - q_fb(k)) - (gam-bet)*m_fb(k);
end

close;
subplot(2,1,1); plot(a,W - W_ns2,a,zeros(1,length(a)),':'); xlabel('Degree of Homogeneity'); ylabel('Diff in Welfare(Open - Close)');
subplot(2,1,2); plot(a,W,a,W_ns2,'--'); 
xlabel('Degree of Homogeneity'); ylabel('Welfare');
legend('Open','Close');


print -djpeg epsFig3

%% Welfare Diff sigma
clear;

bin = 100;

bet = .98; sig = [.01:.01:.5]; thet = .8; phit0 = 1.06; phit1 = 1;
tau = (phit0/phit1) - 1; i_b = .05; g =1.4; a = .5 ; gam = 1 + tau;

for k = 1:length(sig)
    
    q_opt = ((1 + gam/bet)/(g*a))^(1/(a-1));
    q_opt_ns2 = (((gam-bet)/(bet*sig(k)) + 1)/(g*a))^(1/(a-1));
    m_opt = q_opt*sig(k)/(bet*phit1);

    m(k) = q_opt*sig(k)/(bet*phit1);
    m_ns2(k) = q_opt_ns2*bet*phit1;
    
    l(k) = (1-sig(k))*m(k)/sig(k);
    d_b(k) = (((2 + i_b)/(a*g))^(1/(a-1)))/(phit1*bet);
    b_res(k) = max(d_b(k) - m(k),0);
    q_l(k) = bet*phit1*m(k)/sig(k);
    q_b(k) = bet*phit1*(m(k)+b_res(k));
    u = @(q) g*(q^a);
    psi_m(k) = u(q_l(k)) - q_l(k);
    s_b(k) = u(q_b(k)) - q_b(k) - bet*phit1*(1+i_b)*b_res(k);
    x(k) = max((thet/bet*phit1)*(psi_m(k) - bet*phit1*l(k) - s_b(k)),0);
    i(k) = x(k)/l(k);
    
    
    if m(k)/sig(k) < d_b
        b(k) = d_b(k) - m(k)/sig(k);
        q_l_b(k) = bet*phit1*(m(k)/sig(k) + b(k));
        W_b(k) = sig(k)*(u(q_l_b(k)) - q_l_b(k)) - (gam-bet)*m(k);
    else
        b(k) = 0;
        q_l_b(k) = 0;
        W_b(k) = sig(k)*(u(q_b(k)) - q_b(k)) - (gam-bet)*m(k);
    end
    W_l(k) = sig(k)*(u(q_l(k)) - q_l(k)) - (gam-bet)*m(k);
    W(k) = max(W_l(k),W_b(k));
    
    q_ns2(k) = phit1*bet*m_ns2(k);
    W_ns2(k) = (u(q_ns2(k))-q_ns2(k))*sig(k) - (gam-bet)*m_ns2(k);
    
    q_fb(k) = inv(a*g)^(inv(a-1));
    m_fb(k) = q_fb(k)/(bet*phit1);
    W_fb(k) = sig(k)*(u(q_fb(k)) - q_fb(k)) - (gam-bet)*m_fb(k);
end

close;
subplot(2,1,1); plot(sig,W - W_ns2,sig,zeros(1,length(sig)),':'); xlabel('Meeting Prob'); ylabel('Diff in Welfare(Open - Close)');
subplot(2,1,2); plot(sig,W,sig,W_ns2,'--'); 
xlabel('Meeting Prob'); ylabel('Welfare');
legend('Open','Close');


print -djpeg epsFig5
%% 3dplot
clear;

bin = 100;

bet = .98; sig = [.01:.49/100:.5]; thet = .8; phit0 = [1:.2/100:1.2]; phit1 = 1;
tau = (phit0/phit1) - 1; i_b = .05; g = 2; a = [0.01:.49/100:.5] ; gam = 1 + tau;

for j = 1:length(gam)
    for n = 1:length(a)
        for k = 1:length(sig)

            q_opt = ((1 + gam(j)/bet)/(g*a(n)))^(1/(a(n)-1));
            q_opt_ns2 = (((gam(j)-bet)/(bet*sig(k)) + 1)/(g*a(n)))^(1/(a(n)-1));
            m_opt = q_opt*sig(k)/(bet*phit1);

            m(k) = q_opt*sig(k)/(bet*phit1);
            m_ns2(k) = q_opt_ns2*bet*phit1;

            l(k) = (1-sig(k))*m(k)/sig(k);
            d_b(k) = (((2 + i_b)/(a(n)*g))^(1/(a(n)-1)))/(phit1*bet);
            b_res(k) = max(d_b(k) - m(k),0);
            q_l(k) = bet*phit1*m(k)/sig(k);
            q_b(k) = bet*phit1*(m(k)+b_res(k));
            u = @(q) g*q^a(n);
            psi_m(k) = u(q_l(k)) - q_l(k);
            s_b(k) = u(q_b(k)) - q_b(k) - bet*phit1*(1+i_b)*b_res(k);
            x(k) = max((thet/bet*phit1)*(psi_m(k) - bet*phit1*l(k) - s_b(k)),0);
            i(k) = x(k)/l(k);


            if m(k)/sig(k) < d_b
                b(k) = d_b(k) - m(k)/sig(k);
                q_l_b(k) = bet*phit1*(m(k)/sig(k) + b(k));
                W_b(k) = sig(k)*(u(q_l_b(k)) - q_l_b(k)) - (gam(j)-bet)*m(k);
            else
                b(k) = 0;
                q_l_b(k) = 0;
                W_b(k) = sig(k)*(u(q_b(k)) - q_b(k)) - (gam(j)-bet)*m(k);
            end
            W_l(k) = sig(k)*(u(q_l(k)) - q_l(k)) - (gam(j)-bet)*m(k);
            W(k) = max(W_l(k),W_b(k));

            q_ns2(k) = phit1*bet*m_ns2(k);
            W_ns2(k) = (u(q_ns2(k))-q_ns2(k))*sig(k) - (gam(j)-bet)*m_ns2(k);

            q_fb(k) = inv(a(n)*g)^(inv(a(n)-1));
            m_fb(k) = q_fb(k)/(bet*phit1);
            W_fb(k) = sig(k)*(u(q_fb(k)) - q_fb(k)) - (gam(j)-bet)*m_fb(k);
            
            W_store(k,n,j) = W(k)-W_ns2(k);
            if W_store(k,n,j) > 0
                tick(k,n,j) = 1;
                W_ss(k,n,j) = W_store(k,n,j);
            else
                tick(k,n,j) = NaN;
                W_ss(k,n,j) = NaN;
            end
        end
    end
end

close; s = labelvolshow(tick,W_ss)

print -djpeg epsFig6