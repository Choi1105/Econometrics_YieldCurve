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
b21 = 2;

b12 = 3;
b22 = 4;

sig21 = 3;
sig22 = 1;

tau = 123;
p11 = tau-1/tau;

P = zeros(ms, ms);

P(1,1) = p11;
P(1,2) = 1 - p11;
P(2,1) = 0;
P(2,2) = 1;


Ym = zeros(T,1);
xm = 5*rand(T,1);

sm = ones(T,1);
sm(tau+1:T) = 2*ones(T-tau,1);

for t = 1:T
        
    if sm(t) == 1
        b1j = b11;
        b2j = b21;
        sig2j = sig21;
    else
        b1j = b12;
        b2j = b22;
        sig2j = sig22;
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
theta0 = [1;2;3;4;5;6;0.1];

% index
index = [1;2;3;4;5;6;7];
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
