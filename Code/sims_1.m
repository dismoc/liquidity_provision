
%% Bargaining Power

syms x;
format long;
clear;

bin_m = 100; bin_thet = 100;

bet = .99; sig = .25; thet = (0:1/bin_thet:1); phit0 = 1.03; phit1 = 1;
tau = (phit0/phit1) - 1; i_b = .1; g = 1; a = .5; gam = 1 + tau;

q_opt = ((1 + (gam-bet)/(bet*sig))/(g*a))^(1/(a-1));
m_opt = q_opt/(bet*phit1*(1+tau/sig));

q_star = (1/(g*a))^(1/(a-1));
m_star = q_star/(bet*phit1);

k = 1;
j = 1;

eq = (((2 + i_b)/(a*g))^(1/(a-1)))/(phit1*bet);

m = (0:m_opt*sig/(bin_m):m_opt*sig); i = zeros(1,length(m)); l = zeros(1,length(m)); 
b = zeros(1,length(m)); surp = zeros(1,length(m)); check = zeros(1,length(m));

b_hold = zeros(1,length(thet)); l_hold = zeros(1,length(thet)); 
i_hold = zeros(1,length(thet)); m_hold = zeros(1,length(thet)); 
c_hold = zeros(1,length(thet)); surp_hold = zeros(1,length(thet)); 


for k = 1:length(thet)
d_b = (((2 + i_b)/(a*g))^(1/(a-1)))/(phit1*bet);
clear j;
    for j = 1:length(m)
        if (1 + tau/sig + (1-sig)/sig)*m(j) < d_b
            l_max = ((1-sig)/sig)*m(j);
            b_res = max((((2 + i_b)/(a*g))^(1/(a-1)))/(phit1*bet) - (1+tau/sig)*m(j),0);
            
            out_op_max = g*(bet*phit1*((1+tau/sig)*m(j) + l_max))^a -  bet*phit1*((1+tau/sig)*m(j) + l_max) - bet*phit1*(1+i_b)*l_max;
            
            f = @(y) -((g*(bet*phit1*((1+tau/sig)*m(j) + l_max))^a -  bet*phit1*((1+tau/sig)*m(j) + l_max) - bet*phit1*(1+y)*l_max - out_op_max)^(1-thet(k)))*((bet*phit1*l_max*y)^thet(k));
            
            try
            [output] = fmincon(f,[0],[],[],[],[]);
            catch
                l(j) = 0; i(j) = 0;
            end
            i(j) = min(output,i_b); l(j) = l_max;
            
            b(j) = d_b - (1 + tau/sig)*m(j) - l(j);
            surp(j) = g*(phit1*bet*((1 + tau/sig)*m(j) + l(j) + b(j)))^a - phit1*bet*((1 + tau/sig)*m(j) + l(j) + b(j));
            check_2(j) = 0; check(j) = 0;
        else
            % if aggregate money is enough to get optimal q, then no need
            % for public loans.
            l_max = ((1-sig)/sig)*m(j);
            b_res = max((((2 + i_b)/(a*g))^(1/(a-1)))/(phit1*bet) - (1+tau/sig)*m(j),0);
            out_op = g*(bet*phit1*((1+tau/sig)*m(j) + b_res))^a -  bet*phit1*((1+tau/sig)*m(j) + b_res) - bet*phit1*(1+i_b)*b_res;
%           

            % works, first iteration - changes to only include nonlcon
            %f = @(x) -((1+x(2))-(thet(k)/(bet*phit1*x(1)))*(g*(bet*phit1*((1+tau/sig)*m(j) + x(1)))^a -  bet*phit1*((1+tau/sig)*m(j) + x(1)) - out_op));
            
            %Try NBS x1 = l, x2 = i
            f = @(x) -((g*(bet*phit1*((1+tau/sig)*m(j) + x(1)))^a -  bet*phit1*((1+tau/sig)*m(j) + x(1)) - bet*phit1*(1+x(2))*x(1) - out_op)^(1-thet(k)))*((bet*phit1*x(1)*x(2))^thet(k));
            x0 = [l_max, 0.001];
%             nonlcon =@(x) mycon(x,k,phit1,bet,m,j,tau,sig,g,a,out_op,thet);
            try
                [output] = fmincon(f,x0,[],[],[],[],[0 0],[l_max i_b]);
                i(j) = output(2);
                l(j) = output(1);
            catch
                warning('money holdings value higher than private loan.');
                l(j) = 0; i(j) = 0; flag = 1;
            end
            
            d_i = (((2 + i(j))/(a*g))^(1/(a-1)))/(phit1*bet);

            % Check the other constraint in the bilateral problem
            check(j) = ((1-thet(k))*output(1)*phit1*bet*((g*a*(bet*phit1*((1+tau/sig)*m(j) +...
                output(1)))^(a-1))-1-(1+i(j))) + thet(k)*(g*(bet*phit1*((1+tau/sig)*m(j) +...
                output(1)))^a -  bet*phit1*((1+tau/sig)*m(j) + output(1)) - out_op - bet*phit1*(1+i(j))*output(1)));
             
            % check that 1+i = ... first constraint holds 
            check_2(j) = (thet(k)/(bet*phit1*l(j)))*(g*(bet*phit1*((1+tau/sig)*m(j) +...
                l(j)))^a -  bet*phit1*((1+tau/sig)*m(j) + l(j)) - out_op);
            
            
            b(j) = 0;
            
            if (1 + tau/sig)*m(j) > d_i
                surp(j) = g*(phit1*bet*((1 + tau/sig)*m(j) + b(j)))^a - phit1*bet*((1 + tau/sig)*m(j) + b(j));
                l(j) = 0;
            else
                surp(j) = g*(phit1*bet*((1 + tau/sig)*m(j) + l(j) + b(j)))^a - phit1*bet*((1 + tau/sig)*m(j) + l(j) + b(j));
            end
        end
    end
[M, I] = max(l); b_hold(k) = b(I); l_hold(k) = l(I); i_hold(k) = i(I); surp_hold(k) = surp(I)- phit0*m(I); c_hold(k) = check(I); m_hold(k) = m(I);

end

close;
subplot(2,3,1); plot(m,surp-phit0*m,m(I),surp(I)-phit0*m(I),'r.'); ylabel('surplus'); xlabel('Money Holdings');
subplot(2,3,2); plot(m,l,m(I),l(I),'r.'); ylabel('Loan Size'); xlabel('Money Holdings');
subplot(2,3,3); plot(m,b,m(I),b(I),'r.'); ylabel('public loan'); xlabel('Money Holdings');
subplot(2,3,4); plot(m,l.*i,m(I),l(I)*i(I),'r.'); ylabel('Utils to Lender'); xlabel('Money Holdings');
subplot(2,3,5); plot(m,m*(1+tau/sig)+l+b,m(I),m(I)*(1+tau/sig)+l(I)+b(I),'r.'); ylabel('Transfer to Seller'); xlabel('Money Holdings');
subplot(2,3,6); plot(m,i,m(I),i(I),'r.'); ylabel('private interest'); ylim([0 max(i_b)]); xlim([0 max(m)]); xlabel('Money Holdings');

close;
subplot(2,3,[1 2]); plot(thet,surp_hold-phit0*m_hold); ylabel('surplus'); xlabel('Barg Power');
subplot(2,3,5); plot(thet,l_hold); ylabel('private loan'); xlabel('Barg Power');
subplot(2,3,4); plot(thet,i_hold); ylabel('private interest'); xlabel('Barg Power');
subplot(2,3,3); plot(thet,i_hold.*l_hold); ylabel('Utils to Lender'); xlabel('Barg Power');
subplot(2,3,6); plot(thet,m_hold*(1+tau/sig)+l_hold+b_hold); ylabel('d to Seller'); xlabel('Barg Power');


%% Discount Window Rate

syms x;
format long;
clear;

bin_m =100; bin =100;

bet = .99; sig = .25; thet = .75 ; phit0 = 1.02; phit1 = 1;
tau = (phit0/phit1) - 1; i_b = [0:.2/bin:.2]; g = 1; a = .5; gam = 1 + tau;

q_opt = ((2 + (gam-bet)/(bet*sig))/(g*a))^(1/(a-1));
m_opt = q_opt/(bet*phit1*(1+tau/sig + (1-sig)/(sig)));

k = 1;
j = 1;

m = (0:2*m_opt/(bin_m):2*m_opt); i = zeros(1,length(m)); l = zeros(1,length(m)); 
b = zeros(1,length(m)); surp = zeros(1,length(m)); check = zeros(1,length(m));

b_hold = zeros(1,length(i_b)); l_hold = zeros(1,length(i_b)); 
i_hold = zeros(1,length(i_b)); m_hold = zeros(1,length(i_b)); 
c_hold = zeros(1,length(i_b)); surp_hold = zeros(1,length(i_b)); 


for k = 1:length(i_b)
d_b = (((2 + i_b(k))/(a*g))^(1/(a-1)))/(phit1*bet);
clear j;
    for j = 1:length(m)
        if (1 + tau/sig + (1-sig)/sig)*m(j) < d_b
            l_max = ((1-sig)/sig)*m(j);
            b_res = max((((2 + i_b(k))/(a*g))^(1/(a-1)))/(phit1*bet) - (1+tau/sig)*m(j),0);

            out_op_max = g*(bet*phit1*((1+tau/sig)*m(j) + l_max))^a -  bet*phit1*((1+tau/sig)*m(j) + l_max) - bet*phit1*(1+i_b(k))*l_max;
            
            f = @(y) -((g*(bet*phit1*((1+tau/sig)*m(j) + l_max))^a -  bet*phit1*((1+tau/sig)*m(j) + l_max) - bet*phit1*(1+y)*l_max - out_op_max)^(1-thet))*((bet*phit1*l_max*y)^thet);
            
            try
            [output] = fmincon(f,[0],[],[],[],[]);
            i(j) = max(output,i_b(k)); l(j) = l_max;
            catch
                l(j) = 0; i(j) = 0;
            end
                       
            b(j) = d_b - (1 + tau/sig)*m(j) - l(j);
            surp(j) = g*(phit1*bet*((1 + tau/sig)*m(j) + l(j) + b(j)))^a - phit1*bet*((1 + tau/sig)*m(j) + l(j) + b(j));
            check_2(j) = 0; check(j) = 0;
        else
            % if aggregate money is enough to get optimal q, then no need
            % for public loans.
            l_max = ((1-sig)/sig)*m(j);
            b_res = max((((2 + i_b(k))/(a*g))^(1/(a-1)))/(phit1*bet) - (1+tau/sig)*m(j),0);
            out_op = g*(bet*phit1*((1+tau/sig)*m(j) + b_res))^a -  bet*phit1*((1+tau/sig)*m(j) + b_res) - bet*phit1*(1+i_b(k))*b_res;

            %Try NBS x1 = l, x2 = i
            f = @(x) -((g*(bet*phit1*((1+tau/sig)*m(j) + x(1)))^a -  bet*phit1*((1+tau/sig)*m(j) + x(1)) - bet*phit1*(1+x(2))*x(1) - out_op)^(1-thet))*((bet*phit1*x(1)*x(2))^thet);
            
            x0 = [l_max, 0.001];

            try
                [output] = fmincon(f,x0,[],[],[],[],[0 0],[l_max i_b(k)],[],optimset('TolFun', 1e-12));
                i(j) = output(2);
                l(j) = output(1);
                
                check(j) = ((1-thet)*output(1)*phit1*bet*((g*a*(bet*phit1*((1+tau/sig)*m(j) +...
                output(1)))^(a-1))-1-(1+i(j))) + thet*(g*(bet*phit1*((1+tau/sig)*m(j) +...
                output(1)))^a -  bet*phit1*((1+tau/sig)*m(j) + output(1)) - out_op - bet*phit1*(1+i(j))*output(1)));
                l_res =max((((2 + output(2))/(a*g))^(1/(a-1)))/(phit1*bet) - (1+tau/sig)*m(j),0);


            catch
                warning('money holdings value higher than private loan.');
                l(j) = 0; i(j) = 0; flag = 1;
            end
             % Check the other constraint in the bilateral problem
            d_i = (((2 + i(j))/(a*g))^(1/(a-1)))/(phit1*bet);
            liq_prem(j) = g*a*(bet*phit1*((1+tau/sig)*m(j)+l(j)))^(a-1) - 1;
            % check that 1+i = ... first constraint holds 
            
            b(j) = 0;
            
            if (1 + tau/sig)*m(j) > d_i
                surp(j) = g*(phit1*bet*((1 + tau/sig)*m(j) + b(j)))^a - phit1*bet*((1 + tau/sig)*m(j) + b(j));
                l(j) = 0;
            else
                surp(j) = g*(phit1*bet*((1 + tau/sig)*m(j) + l(j) + b(j)))^a - phit1*bet*((1 + tau/sig)*m(j) + l(j) + b(j));
            end
        end
        
        buyer_surp(j) = sig*(surp(j)/phit1 - gam*m(j) - bet*(1+i(j))*l(j));
        lender_surp(j) = (bet*(i(j)*l(j)+m(j)) - gam*m(j))*(1-sig);
    end
W = lender_surp + buyer_surp;
[M, I] = max(l); b_hold(k) = b(I); l_hold(k) = l(I); i_hold(k) = i(I); surp_hold(k) = surp(I); c_hold(k) = check(I); m_hold(k) = m(I);

end

close;
subplot(2,3,1); plot(m,W,m(I),W(I),'r.'); ylabel('surplus'); xlabel('Money Holdings');
subplot(2,3,2); plot(m,l,m(I),l(I),'r.'); ylabel('Loan Size'); xlabel('Money Holdings');
subplot(2,3,3); plot(m,b,m(I),b(I),'r.'); ylabel('public loan'); xlabel('Money Holdings');
subplot(2,3,4); plot(m,l.*i,m(I),l(I)*i(I),'r.'); ylabel('Utils to Lender'); xlabel('Money Holdings');
subplot(2,3,5); plot(m,m*(1+tau/sig)+l+b,m(I),m(I)*(1+tau/sig)+l(I)+b(I),'r.'); ylabel('Transfer to Seller'); xlabel('Money Holdings');
subplot(2,3,6); plot(m,i,m(I),i(I),'r.'); ylabel('private interest'); ylim([0 max(i_b)]); xlim([0 max(m)]); xlabel('Money Holdings');


W = surp_hold-phit0*m_hold;
close;
subplot(2,3,1); plot(i_b,W); ylabel('Surplus'); xlabel('DW Rate');
subplot(2,3,2); plot(i_b,m_hold); ylabel('Money Holdings'); xlabel('DW Rate');
subplot(2,3,5); plot(i_b,l_hold); ylabel('private loan'); xlabel('DW Rate');
subplot(2,3,4); plot(i_b,i_hold); ylabel('private interest'); xlabel('DW Rate');
subplot(2,3,3); plot(i_b,i_hold.*l_hold); ylabel('Utils to Lender'); xlabel('DW Rate');
subplot(2,3,6); plot(i_b,b_hold); ylabel('Public Borrowing'); xlabel('DW Rate');

clear i l; i = [0:i_b(end)/100:i_b(end)]; l=[0:l_max/100:l_max];
for k = 1:101
    for n=1:101
        z(n,k) = f([i(n) l(k)]);
    end
end
hold on; [M I] = max(z)
surf(i,l,z); surf(
for n=1:360
    camroll(1)
    drawnow
end

%% Inflation
syms x;
format long;
clear;

bin_m = 100; bin = 50;

bet = .99; sig = .25; thet = .75 ; phit1 = 1;
i_b = .1; g = 1; a = .5; gam = [1:.2/bin:1.2];


k = 1;
j = 1;

i = zeros(1,bin_m); l = zeros(1,bin_m); 
b = zeros(1,bin_m); surp = zeros(1,bin_m); check = zeros(1,bin_m);

b_hold = zeros(1,length(i_b)); l_hold = zeros(1,length(i_b)); 
i_hold = zeros(1,length(i_b)); m_hold = zeros(1,length(i_b)); 
c_hold = zeros(1,length(i_b)); surp_hold = zeros(1,length(i_b)); 


for k = 1:length(gam)
d_b = (((2 + i_b)/(a*g))^(1/(a-1)))/(phit1*bet);
tau = gam(k) - 1;
q_opt = ((2 + (gam(k)-bet)/(bet*sig))/(g*a))^(1/(a-1));
m_opt = q_opt/(bet*phit1*(1+tau/sig + (1-sig)/(sig)));
m = (0:2*m_opt/(bin_m):2*m_opt);
phit0 = gam;

clear j;
    for j = 1:length(m)
        if (1 + tau/sig + (1-sig)/sig)*m(j) < d_b
            l_max = ((1-sig)/sig)*m(j);
            b_res = max((((2 + i_b)/(a*g))^(1/(a-1)))/(phit1*bet) - (1+tau/sig)*m(j),0);

            out_op_max = g*(bet*phit1*((1+tau/sig)*m(j) + l_max))^a -  bet*phit1*((1+tau/sig)*m(j) + l_max) - bet*phit1*(1+i_b)*l_max;
            
            f = @(y) -((g*(bet*phit1*((1+tau/sig)*m(j) + l_max))^a -  bet*phit1*((1+tau/sig)*m(j) + l_max) - bet*phit1*(1+y)*l_max - out_op_max)^(1-thet))*((bet*phit1*l_max*y)^thet);
            
            try
            [output] = fmincon(f,[0],[],[],[],[]);
            catch
                l(j) = 0; i(j) = 0;
            end
            i(j) = output(end); l(j) = l_max;            
            b(j) = d_b - (1 + tau/sig)*m(j) - l(j);
            surp(j) = g*(phit1*bet*((1 + tau/sig)*m(j) + l(j) + b(j)))^a - phit1*bet*((1 + tau/sig)*m(j) + l(j) + b(j));
            check_2(j) = 0; check(j) = 0;
        else
            % if aggregate money is enough to get optimal q, then no need
            % for public loans.
            l_max = ((1-sig)/sig)*m(j);
            b_res = max((((2 + i_b)/(a*g))^(1/(a-1)))/(phit1*bet) - (1+tau/sig)*m(j),0);
            out_op = g*(bet*phit1*((1+tau/sig)*m(j) + b_res))^a -  bet*phit1*((1+tau/sig)*m(j) + b_res) - bet*phit1*(1+i_b)*b_res;
%           

            % works, first iteration - changes to only include nonlcon
            %f = @(x) -((1+x(2))-(thet(k)/(bet*phit1*x(1)))*(g*(bet*phit1*((1+tau/sig)*m(j) + x(1)))^a -  bet*phit1*((1+tau/sig)*m(j) + x(1)) - out_op));
            
            %Try NBS x1 = l, x2 = i
            f = @(x) -((g*(bet*phit1*((1+tau/sig)*m(j) + x(1)))^a -  bet*phit1*((1+tau/sig)*m(j) + x(1)) - bet*phit1*(1+x(2))*x(1) - out_op)^(1-thet))*((bet*phit1*x(1)*x(2))^thet);
            
            x0 = [l_max, 0];
%             nonlcon =@(x) mycon(x,k,phit1,bet,m,j,tau,sig,g,a,out_op,thet);
            try
                [output] = fmincon(f,x0,[],[],[],[],[0 0],[l_max i_b]);
                i(j) = output(2);
                l(j) = output(1);
            catch
                warning('money holdings value higher than private loan.');
                l(j) = 0; i(j) = 0; flag = 1;
            end
%             
             d_i = (((2 + i(j))/(a*g))^(1/(a-1)))/(phit1*bet);

             % Check the other constraint in the bilateral problem
            check(j) = ((1-thet)*output(1)*phit1*bet*((g*a*(bet*phit1*((1+tau/sig)*m(j) +...
                output(1)))^(a-1))-1-(1+i(j))) + thet*(g*(bet*phit1*((1+tau/sig)*m(j) +...
                output(1)))^a -  bet*phit1*((1+tau/sig)*m(j) + output(1)) - out_op - bet*phit1*(1+i(j))*output(1)));
            %without nonlcon
             %l(j) = max(min(l_max,l_res),0) ;
            %with nonlcon
             %l(j) = max(min(l_max,output(1)),0);
             
            % check that 1+i = ... first constraint holds 
            check_2(j) = (thet/(bet*phit1*l(j)))*(g*(bet*phit1*((1+tau/sig)*m(j) +...
                l(j)))^a -  bet*phit1*((1+tau/sig)*m(j) + l(j)) - out_op);
            
            
            b(j) = 0;
            
            if (1 + tau/sig)*m(j) > d_i
                surp(j) = g*(phit1*bet*((1 + tau/sig)*m(j) + b(j)))^a - phit1*bet*((1 + tau/sig)*m(j) + b(j));
                l(j) = 0;
            else
                surp(j) = g*(phit1*bet*((1 + tau/sig)*m(j) + l(j) + b(j)))^a - phit1*bet*((1 + tau/sig)*m(j) + l(j) + b(j));
            end
            
            if l(j) >= l_max
                pl_con(j) = 1;
            else
                pl_con(j) = 0;
            end
        end
    end
    
[M, I] = max(l); b_hold(k) = b(I); l_hold(k) = l(I); i_hold(k) = i(I); surp_hold(k) = surp(I); c_hold(k) = check(I); m_hold(k) = m(I); 

pri_loan_use(k) = pl_con(I);

end

close;
subplot(2,3,1); plot(m,surp-gam(k)*m,m(I),surp(I)-gam(k)*m(I),'r.'); ylabel('surplus'); xlabel('Money Holdings');
subplot(2,3,2); plot(m,l,m(I),l(I),'r.'); ylabel('Loan Size'); xlabel('Money Holdings');
subplot(2,3,3); plot(m,b,m(I),b(I),'r.'); ylabel('public loan'); xlabel('Money Holdings');
subplot(2,3,4); plot(m,l.*i,m(I),l(I)*i(I),'r.'); ylabel('Utils to Lender'); xlabel('Money Holdings');
subplot(2,3,5); plot(m,m*(1+tau/sig)+l+b,m(I),m(I)*(1+tau/sig)+l(I)+b(I),'r.'); ylabel('Transfer to Seller'); xlabel('Money Holdings');
subplot(2,3,6); plot(m,i,m(I),i(I),'r.'); ylabel('private interest'); ylim([0 max(i_b)]); xlim([0 max(m)]); xlabel('Money Holdings');

close;
subplot(2,3,1); plot(gam-1,surp_hold - gam.*m_hold); ylabel('surplus'); xlabel('Inflation');
subplot(2,3,2); plot(gam-1,m_hold); ylabel('Money Holdings'); xlabel('Inflation');
subplot(2,3,5); plot(gam-1,l_hold); ylabel('private loan'); xlabel('Inflation');
subplot(2,3,4); plot(gam-1,i_hold); ylabel('private interest'); xlabel('Inflation');
subplot(2,3,3); plot(gam-1,i_hold.*l_hold); ylabel('Utils to Lender'); xlabel('Inflation');
subplot(2,3,6); plot(gam-1,m_hold*(1+tau/sig)+l_hold+b_hold); ylabel('d to Seller'); xlabel('Inflation');

