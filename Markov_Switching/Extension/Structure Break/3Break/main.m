%% Structure Break
% y(t) = x(t)*B(St) + e(t)
% e(t) ~ N(0,sig2(St))
% B(St) = B0(1-St)+B1(St)
% Sig2
% St = 0 or 1
clear;
clc;

%% Step 1: DGP %%
T = 1002;
ms = 3;

b11 = 1;
b12 = 2;
b21 = 3;
b22 = 4;
b13 = 5;
b23 = 6;

sig21 = 3;
sig22 = 2;
sig23 = 1;

tau1 = 123;
tau2 = 456;

p11 = (tau1 - 1)/tau1;
p22 = (tau2 - tau1 - 1)/(tau2 - tau1);

P = zeros(ms, ms);
P(1,1) = p11;
P(1,2) = 1 - p11;
P(1,3) = 0; % jump 불가능.
P(2,1) = 0;
P(2,2) = p22;
P(2,3) = 1 - p22;
P(3,1) = 0;
P(3,2) = 0;
P(3,3) = 1; % 뒤로 불가능.


Ym = zeros(T,1);
xm = 5*rand(T,1);

sm = ones(T,1);
sm(tau1+1:tau2) = 2*ones(tau2-tau1,1);
sm(tau2+1:T) = 3*ones(T-tau2,1);

for t = 1:T
        
   if sm(t) == 1
        b1j = b11;
        b2j = b21;
        sig2j = sig21;
    elseif sm(t) == 2
        b1j = b12;
        b2j = b22;
        sig2j = sig22;
    else
        b1j = b13;
        b2j = b23;
        sig2j = sig23;
    end

    % Generate y(t)
    Ym(t) = b1j + b2j*xm(t) + sqrt(sig2j)*randn(1,1);

end
%% Step 2: Estimation %%
% Data
Y = Ym; 
X = [ones(T,1) xm]; 
Data = [Y X];
T = rows(X);
k = cols(X);

% initial value 
theta0 = [1;2;3;4;5;6;1;2;3;0.5;0.5];

% index
index = [1;2;3;4;5;6;7;8;9;10;11];
printi = 1;

% Optimization
[thetamx, fmax, Vj, Vinv] = SA_Newton(@lnlik,@paramconst,theta0,Data,printi,index);
[~, Spm] = lnlik(thetamx, Data);



%% Step 3: Results %%
i = 1:T;
i = i';

smm = sm(1:end);

plot(i, abs(smm-2),  'LineStyle', '-', 'LineWidth', 2.0, 'Color', 'b');
hold on
plot(i, Spm(:,1),  'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'r');
plot(i, 1.01*ones(T,1),  'LineStyle', ':', 'LineWidth', 1.0, 'Color', 'k');
plot(i, -0.01*ones(T,1),  'LineStyle', ':', 'LineWidth', 1.0, 'Color', 'k');
legend('DGP', 'Estimates')
title('Prob. of Regime 1');
