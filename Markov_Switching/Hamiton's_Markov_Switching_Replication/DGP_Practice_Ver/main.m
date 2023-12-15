%% Application 1 : Hamilton's Markov-Switching Model of Business Fluctuations
% Textbook : 78p ~ 81p

% Y(t) = log(y_t) - log(y_t-1) 
% Y(t) - Mu(St) = phi1*(Y(t-1) - Mu(St-1)) + phi2*(Y(t-2) - Mu(St-2))  
% + phi3*(Y(t-3) - Mu(St-3)) + phi4*(Y(t-4) - Mu(St-4)) +  e(t)

% e(t) ~ iidN(0, Sig2)
% Mu(St) = Mu0(St-1) + Mu1(St)

% Pr[St = 1 | St-1 = 1] = p
% Pr[St = 1 | St-1 = 0] = 1 - p
% Pr[St = 0 | St-1 = 0] = q
% Pr[St = 0 | St-1 = 1] = 1 - q

clear;
clc;

%% Step 1: DGP %%
T = 1504;
Regime = 2;

% Parameter
b11 = 1;
b12 = 6;

phi1 = 0.5;
phi2 = -0.1;
phi3 = -0.3;
phi4 = 0.2;

sig2 = 0.1;

p0 = 0.95;
q0 = 0.95;

% Prob-Matrix
P = zeros(Regime, Regime);

P(1,1) = p0;
P(1,2) = 1 - p0;
P(2,1) = 1 - q0;
P(2,2) = q0;

Im = ones(Regime,1);
Zm = zeros(Regime,1);

A = [eye(Regime) - P' ; Im'];
SS_prob = inv(A'*A)*A'*[Zm;1];

pr1 = SS_prob(1);
pr2 = SS_prob(2);

% Pre-allocation
sm = zeros(T, 1);
Ym = zeros(T, 1);
Lb1j = zeros(T, 1);

% Ym(1) Generate Code (Assumption : St = 1)
sm(1) = 1;
Lb1j(1) = b11;
Ym(1) = b11 + sqrt(sig2)*randn(1,1);

% Ym(2) Generate Code (Assumption : St = 1)
sm(2) = 1;
Lb1j(2) = b11;
Ym(2) = b11 + phi1*(Ym(1)-b11) + sqrt(sig2)*randn(1,1);

% Ym(3) Generate Code (Assumption : St = 1)
sm(3) = 1;
Lb1j(3) = b11;
Ym(3) = b11 + phi1*(Ym(2)-b11) + phi2*(Ym(1)-b11) + sqrt(sig2)*randn(1,1);

% Ym(4) Generate Code (Assumption : St = 1)
sm(4) = 1;
Lb1j(4) = b11;
Ym(4) = b11 + phi1*(Ym(3)-b11) + phi2*(Ym(2)-b11) + phi3*(Ym(1)-b11) + sqrt(sig2)*randn(1,1);

% Ym(t) Generate Code (Assumpthon : St(5) = 1)
Lst = 1;

for t = 5:T
    if Lst == 1
        pr1 = P(1,1);
    else
        pr1 = P(2,1);
    end

    if rand(1, 1) < pr1
        st = 1;
    else
        st = 2;
    end

    sm(t, 1) = st;

    if st == 1
        b1j = b11;
    else
        b1j = b12;
    end

    % Generate y(t)
    Ym(t) = b1j + phi1*(Ym(t-1)-Lb1j(t-1,1)) + phi2*(Ym(t-2)-Lb1j(t-2,1)) + phi3*(Ym(t-3)-Lb1j(t-3,1)) + phi4*(Ym(t-4)-Lb1j(t-4,1)) + sqrt(sig2)*randn(1,1);

    % Updating
    Lst = st;
    Lb1j(t,1) = b1j;
end

% Burn-in
Ym_New = Ym(501:end,1);
sm_New = sm(501:end,1);
nT = rows(Ym_New);

%% Step 2: Estimation %%
% Data
Y = Ym_New(5:end); 
X = [Ym_New(1:nT-4) Ym_New(2:nT-3) Ym_New(3:nT-2) Ym_New(4:nT-1)]; 
Data = [Y X];
T = rows(X);
k = cols(X);

% initial value 
% Beta1 Beta2 Phi1 Phi2 Phi3 Phi4 Sig2 p0 q0
theta0 = [1.2 ; 3 ; 0.5 ;-0.1; -0.7 ; 0.2 ; 1 ; 0.95 ; 0.95];

% index
index = [1;2;3;4;5;6;7;8;9];
printi = 1;

% Optimization
[thetamx, fmax, Vj, Vinv] = SA_Newton(@lnlik,@paramconst,theta0,Data,printi,index);
[~, Spm] = lnlik(thetamx, Data);


%% Step 3: Results %%
i = 1:T;
i = i';

smm = sm_New(5:end);

plot(i, abs(smm-2),  'LineStyle', '-', 'LineWidth', 2.0, 'Color', 'b');
hold on
plot(i, Spm(:,1),  'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'r');
plot(i, 1.01*ones(T,1),  'LineStyle', ':', 'LineWidth', 1.0, 'Color', 'k');
plot(i, -0.01*ones(T,1),  'LineStyle', ':', 'LineWidth', 1.0, 'Color', 'k');
legend('DGP', 'Estimates')
title('Prob. of Regime 1');
