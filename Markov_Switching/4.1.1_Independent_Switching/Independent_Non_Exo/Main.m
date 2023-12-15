%% Maximum likelihood estimation
% Independent Switching
% y(t) = x(t)*B(St) + e(t)
% e(t) ~ N(0,sig2(St))
% B(St) = B0(1-St)+B1(St)
% Sig2
% St = 0 or 1
clear;
clc;

%% Step 1: DGP %%
T = 2000;

% p0 값 설정
p0 = 2;
Proba = (exp(p0))/(1+exp(p0));
% St 결정
St = zeros(1, T); % 초기값 0으로 설정
for i = 1:T
    if rand < exp(p0)/(1+exp(p0))
        St(i) = 1;
    end
end
St = St';

% Sig2
Sig2_1 = 2;
Sig2_2 = 1;

% Beta_0 and Beta_1
Beta_1 = 15;
Beta_2 = 1;

Ym = zeros(T,1);
X1m = 2*rand(T,1);
em = zeros(T,1);
emx = randn(T,1);

for t = 1:T
    em(t) = sqrt(Sig2_1)*emx(t)*(1-St(t)) + sqrt(Sig2_2)*emx(t)*St(t);
end

for t = 1:T
    Ym(t) = X1m(t)*Beta_1*(1-St(t)) + X1m(t)*Beta_2*St(t) + em(t);
end

%% Step 2: Estimation %%
% Data
Y = Ym; 
X = X1m; 
Data = [Y X St];
T = rows(X);
k = cols(X);

% initial value 
theta0 = [15 ; 2 ; 1 ; 2 ; 0.1 ];

% index
index = [1;2;3;4;5];
printi = 1;

% Optimization
[thetamx, fmax, Vj, Vinv] = SA_Newton(@lnlik,@paramconst,theta0,Data,printi,index);

% Compute t-value and p-value
diag_cov = diag(Vj);
stde = sqrt(diag_cov);
t_val = thetamx./stde; 
p_val = 2*(1 - cdf('t', abs(t_val), T-k)); 

%% Step 3: Results %%
disp('=========================================');
disp(['  Estimates ','  t value ', '  p value ']); 
disp('=========================================');
disp([thetamx t_val p_val]); 
disp('=========================================');
disp(["p0로 계산된 p 출력",exp(p0)/(1+exp(p0))])
disp(["St=1인 비율 출력",sum(St)/T]); 
          