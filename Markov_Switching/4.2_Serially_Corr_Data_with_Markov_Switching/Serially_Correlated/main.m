%% Serially Correlated and Markov Switching
% Textbook : 64p ~ 68p

% Y(t) - Mu(St) = phi1*(Y(t-1) - Mu(St-1)) + e(t), e(t) ~ N(0, Sig2(St))
clear;
clc;

%% Step 1: DGP %%
T = 1502;
ms = 2;

b11 = 3;
b12 = 5;

phi1 = 0.9;

sig21 = 1;
sig22 = 5;

p0 = 0.95;
q0 = 0.95;

P = zeros(ms, ms);

P(1,1) = p0;
P(1,2) = 1 - p0;
P(2,1) = 1 - q0;
P(2,2) = q0;

pr1 = (1 - p0) / (2 - p0 - q0);
pr2 = (1 - q0) / (2 - p0 - q0);

Lst = 2; % 초기값 설정

sm = zeros(T, 1);
Ym = zeros(T, 1);
Lb1j = zeros(T, 1);
Lb1j(1, 1) = b11;

Ym(1) = phi1*b12 + sqrt(sig22)*randn(1, 1);

for t = 2:T
    if Lst == 1
        pr1 = p0;
    else
        pr1 = 1 - q0;
    end

    if rand(1, 1) < pr1
        st = 1;
    else
        st = 2;
    end

    sm(t, 1) = st;

    if st == 1
        b1j = b11;
        sig2j = sig21;
    else
        b1j = b12;
        sig2j = sig22;
    end

    % Generate y(t)
    Ym(t) = phi1*(Ym(t-1) - Lb1j(t-1)) + b1j + sqrt(sig2j) * randn(1, 1);

    % Updating
    Lst = st;
    Lb1j(t,1) = b1j;
end

Y0 = Ym(1002,1);
Ym = Ym(503:1502,1);
sm = sm(503:1502,1);

T = rows(Ym);
%% Step 2: Estimation %%
% Data
Y = Ym(2:end); 
YX = [Ym(1:T-1)];
X = [ones(999,1) YX]; 
Data = [Y X];
T = rows(X);
k = cols(X);

% initial value 
theta0 = [1 ; 5 ; 0.5 ;0.1; 0.4 ; 0.9 ; 0.9];

% index
index = [1;2;3;4;5;6;7];
printi = 1;

% Optimization
[thetamx, fmax, Vj, Vinv] = SA_Newton(@lnlik,@paramconst,theta0,Data,printi,index);
[~, Spm] = lnlik(thetamx, Data);



%% Step 3: Results %%
i = 1:T;
i = i';

smm = sm(2:end);

plot(i, abs(smm-2),  'LineStyle', '-', 'LineWidth', 2.0, 'Color', 'b');
hold on
plot(i, Spm(:,1),  'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'r');
plot(i, 1.01*ones(T,1),  'LineStyle', ':', 'LineWidth', 1.0, 'Color', 'k');
plot(i, -0.01*ones(T,1),  'LineStyle', ':', 'LineWidth', 1.0, 'Color', 'k');
legend('DGP', 'Estimates')
title('Prob. of Regime 1');
