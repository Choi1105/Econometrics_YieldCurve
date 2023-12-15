%% Markov Switching with different regime

clear;
clc;

%% Step 1: DGP %%

T = 1500;
RG1 = 2;
RG2 = 2;

b11 = 5;
b12 = 0.1;

sig21 = 5;
sig22 = 0.1;

% RG1 Prob
p1 = 0.97;
q1 = 0.97;

P1 = zeros(RG1, RG1);

P(1,1) = p1;
P(1,2) = 1 - p1;
P(2,1) = 1 - q1;
P(2,2) = q1;

pr11 = (1 - p1) / (2 - p1 - q1);
pr12 = (1 - q1) / (2 - p1 - q1);

% RG2 Prob
p2 = 0.95;
q2 = 0.95;

P2 = zeros(RG2, RG2);

P1(1,1) = p2;
P1(1,2) = 1 - p2;
P1(2,1) = 1 - q2;
P1(2,2) = q2;

pr21 = (1 - p1) / (2 - p1 - q1);
pr22 = (1 - q1) / (2 - p1 - q1);

qm = zeros(T,1);
sm = zeros(T,1);
Ym = zeros(T,1);

Lst = 1;
Qst = 2;

Ym(1) = b11 + sqrt(sig22)*randn(1,1);


for t = 2:T

    if Lst == 1   
        pr1 = P(1,1);
    else          
        pr1 = P(2,1);
    end

    if Qst == 1
        pr11 = P1(1,1);
    else
        pr11 = P1(2,1);
    end

    if rand(1,1) < pr1
        st = 1;
    else                 
        st = 2;
    end

    if rand(1,1) < pr11
        qt = 1;
    else
        qt = 2;
    end

qm(t,1) = qt;
sm(t,1) = st;
        
    if st == 1
        b1j = b11;
    else
        b1j = b12;
    end
    
    if qt == 1
        sig2j = sig21;
    else
        sig2j = sig22;
    end

    % Generate y(t)
    Ym(t) = b1j + sqrt(sig2j)*randn(1,1);

    % Updating
    Lst = st;
    Qst = qt;
end

Ym = Ym(501:end,1);

%% Step 2: Estimation %%
% Data
Y = Ym; 
X = [ones(1000,1)]; 
Data = [Y X];
T = rows(X);
k = cols(X);

% initial value 
theta0 = [5 ; 0.1 ; 5 ; 0.1 ; 0.97 ; 0.97 ; 0.97 ; 0.97];

% index
index = [1;2;3;4;5;6;7;8];
printi = 1;

% Optimization
[thetamx, fmax, Vj, Vinv] = SA_Newton(@lnlik,@paramconst,theta0,Data,printi,index);
[~, Spm, Rpm] = lnlik(thetamx, Data);

%% Step 3: Results %%
i = 1:T;
i = i';

sm = sm(501:end,1);
qm = qm(501:end,1);

figure
plot(i, abs(sm-2),  'LineStyle', '-', 'LineWidth', 2.0, 'Color', 'b');
hold on
plot(i, Spm(:,1),  'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'r');
plot(i, 1.01*ones(T,1),  'LineStyle', ':', 'LineWidth', 1.0, 'Color', 'k');
plot(i, -0.01*ones(T,1),  'LineStyle', ':', 'LineWidth', 1.0, 'Color', 'k');
legend('DGP', 'Estimates')
title('Prob. of Regime S1');

figure
plot(i, abs(qm-2),  'LineStyle', '-', 'LineWidth', 2.0, 'Color', 'b');
hold on
plot(i, Rpm(:,1),  'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'r');
plot(i, 1.01*ones(T,1),  'LineStyle', ':', 'LineWidth', 1.0, 'Color', 'k');
plot(i, -0.01*ones(T,1),  'LineStyle', ':', 'LineWidth', 1.0, 'Color', 'k');
legend('DGP', 'Estimates')
title('Prob. of Regime R1');