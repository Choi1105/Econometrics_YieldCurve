%% Markov Switching with 3 Regimes %%

% Y(t) = beta1 + beta2*x(t) + e, e ~ N(0 , sig^2)

clear;
clc;

%% Step 1: DGP %%
T = 1002;
ms = 3;

b11 = 1;
b12 = 2;
b13 = 3;
b21 = 1;
b22 = 2;
b23 = 3;

sig21 = 0.1;
sig22 = 2;
sig23 = 5;

p11 = 0.9;
p12 = 0.05; % else = 0.3
p21 = 0.03;
p22 = 0.9; % else = 0.4
p31 = 0.03;
p32 = 0.02; % else = 0.4
% 위 값이 너무 낮으면 DGP 제대로 나오지 않음.

P = zeros(ms, ms);

P(1,1) = p11;
P(1,2) = p12;
P(1,3) = 1 - p11 - p12;
P(2,1) = p21;
P(2,2) = p22;
P(2,3) = 1 - p21 - p22;
P(3,1) = p31;
P(3,2) = p32;
P(3,3) = 1 - p31 - p32;

Lst = 1; % 초기값 설정

sm = zeros(T, 1);
Ym = zeros(T, 1);
xm = 5*rand(T,1);
Lb1j = zeros(T, 1);
Lb1j(1, 1) = b11;

Ym(1) = b11 + b21*1 + sqrt(sig21)*randn(1,1);

% Pr지정 시작
for t = 2:T
    if Lst == 1
        pr1 = P(1,1);
    elseif Lst ==2
        pr1 = P(2,1);
    else 
        pr1 = P(3,1);
    end

    if Lst == 1
        pr2 = P(1,2);
    elseif Lst ==2
        pr2 = P(2,2);
    else 
        pr2 = P(3,2);
    end

% Pr 지정 끝
% St 지정 -> 경우 3가지
    random = rand(1,1); % 1 iter안에서 바뀌면 안됨.

    if random < pr1
        st = 1;
    elseif pr1 <= random && random < pr1 + pr2
        st = 2;
    else
        st = 3;
    end

    sm(t, 1) = st;

    if st == 1
        b1j = b11;
        b2j = b21;
        sig2j = sig21;
    elseif st == 2
        b1j = b12;
        b2j = b22;
        sig2j = sig22;
    else
        b1j = b13;
        b2j = b23;
        sig2j = sig23;
    end

    % Generate y(t)
    Ym(t) =b1j + b2j*xm(t) + sqrt(sig2j)*randn(1,1);

    % Updating
    Lst = st;
    
end
%Y0 = Ym(502,1);
%Ym = Ym(503:1002,1);
%sm = sm(503:1002,1);
%xm = xm(503:1002,1);

T = rows(Ym);
%% Step 2: Estimation %%
% Data
Y = Ym; 
X = [ones(1002,1) xm]; 
Data = [Y X];
T = rows(X);
k = cols(X);

% initial value 
theta0 = [b11;b21;b12;b22;b13;b23;sig21;sig22;sig23;p11;p12;p21;p22;p31;p32];

% index
index = [1;2;3;4;5;6;7;8;9;10;11;12;13;14;15];
printi = 1;

% Optimization
[thetamx, fmax, Vj, Vinv] = SA_Newton(@lnlik,@paramconst,theta0,Data,printi,index);
[~, Spm] = lnlik(thetamx, Data);



%% Step 3: Results %%
i = 1:T;
i = i';

smm1 = zeros(T,1);
smm2 = zeros(T,1);
smm3 = zeros(T,1);

for t = 1:T
    if sm(t) == 1
        smm1(t) = 1;
    elseif sm(t) == 2
        smm2(t) = 1;
    else
        smm3(t) = 1;
    end
end

%% Step 4: Display %%
figure
plot(i, abs(smm1),  'LineStyle', '-', 'LineWidth', 2.0, 'Color', 'b');
hold on
plot(i, Spm(:,1),  'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'r');
plot(i, 1.01*ones(T,1),  'LineStyle', ':', 'LineWidth', 1.0, 'Color', 'k');
plot(i, -0.01*ones(T,1),  'LineStyle', ':', 'LineWidth', 1.0, 'Color', 'k');
legend('DGP', 'Estimates')
title('Prob. of Regime 1');

figure
plot(i, abs(smm2),  'LineStyle', '-', 'LineWidth', 2.0, 'Color', 'b');
hold on
plot(i, Spm(:,2),  'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'r');
plot(i, 1.01*ones(T,1),  'LineStyle', ':', 'LineWidth', 1.0, 'Color', 'k');
plot(i, -0.01*ones(T,1),  'LineStyle', ':', 'LineWidth', 1.0, 'Color', 'k');
legend('DGP', 'Estimates')
title('Prob. of Regime 1');

figure
plot(i, abs(smm3),  'LineStyle', '-', 'LineWidth', 2.0, 'Color', 'b');
hold on
plot(i, Spm(:,3),  'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'r');
plot(i, 1.01*ones(T,1),  'LineStyle', ':', 'LineWidth', 1.0, 'Color', 'k');
plot(i, -0.01*ones(T,1),  'LineStyle', ':', 'LineWidth', 1.0, 'Color', 'k');
legend('DGP', 'Estimates')
title('Prob. of Regime 1');



