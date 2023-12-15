%% OLS estimation
% Non Linear Model
% y(t) = x(t)*B(St) + e(t)
% e(t) ~ N(0,sig2(St))
% B(St) = B0(1-St)+B1(St)
% Sig2
% St = 0 or 1
clear;
clc;
%% Step 1: DGP %%
T = 10000;

St = zeros(T,1);
for v = [1 5 8]
    St(v*1000:(v*1000)+1000,1) = 1;
end

% Sig2
Sig2_1 = 0.1;
Sig2_2 = 0.5;

% Beta_0 and Beta_1
Beta_0 = 1.2;
Beta_1 = 5.9;

Ym = zeros(T,1);
X1m = 5*rand(T,1);
em = zeros(T,1);
emx = randn(T,1);

for t = 1:10000
    em(t) = sqrt(Sig2_1)*emx(t)*(1-St(t)) + sqrt(Sig2_2)*emx(t)*St(t);
end

for t = 1:10000
    Ym(t) = X1m(t)*Beta_0*(1-St(t)) + X1m(t)*Beta_1*St(t) + em(t);
end

%% Step 2: Estimation %%
k = 1;

Ymm1 = Ym([1:999,2001:4999,6001:7999,9001:10000],1);
Ymm2 = Ym([1000:2000,5000:6000,8000:9000],1);

T1 = rows(Ymm1);
T2 = rows(Ymm2);

X1 = X1m([1:999,2001:4999,6001:7999,9001:10000],1);
X2 = X1m([1000:2000,5000:6000,8000:9000],1);

bhat1 = inv(X1'*X1)*X1'*Ymm1;
bhat2 = inv(X2'*X2)*X2'*Ymm2;

Yhat1 = X1*bhat1;
Yhat2 = X2*bhat2;

ehat1 = Ymm1 - Yhat1;
ehat2 = Ymm2 - Yhat2;

sig2hat1 = ehat1'*ehat1/(T1-k);
sig2hat2 = ehat2'*ehat2/(T2-k);

varbhat1 = sig2hat1*inv(X1'*X1);
varbhat2 = sig2hat2*inv(X2'*X2);

stde1 = sqrt(diag(varbhat1));
stde2 = sqrt(diag(varbhat2));

B0 = zeros(k,1);

t_val1 = (bhat1 - B0)./stde1;
t_val2 = (bhat2 - B0)./stde2;

p_val1 = 2*(1 - cdf("t",abs(t_val1),T1-k));
p_val2 = 2*(1 - cdf("t",abs(t_val2),T2-k));


RSS1 = ehat1'*ehat1;
RSS2 = ehat2'*ehat2;

TSS1 = (Ymm1 - mean(Ymm1))'*(Ymm1 - mean(Ymm1));
TSS2 = (Ymm2 - mean(Ymm2))'*(Ymm2 - mean(Ymm2));

R21 = 1- RSS1/TSS1;
R22 = 1- RSS2/TSS2;

R2_1 = 1 - (RSS1*(T1-1))/(TSS1*(T1-k));
R2_2 = 1 - (RSS2*(T2-1))/(TSS2*(T2-k));

SC1 = log(RSS1/T1) - k*log(T1)/T1;
SC2 = log(RSS2/T2) - k*log(T2)/T2;

AIC1 = log(RSS1/T1) - 2*k/T1;
AIC2 = log(RSS2/T2) - 2*k/T2;

%%
disp("==============================================================================================")
disp(["       St = 1 OLS"])
disp("        Bhat1           Sig2hat1          varbhat          stde                 tval                  pval          R^2          R^2bar           SC            AIC")
disp([bhat1 sig2hat1 varbhat1 stde1 t_val1 p_val1 R21 R2_1 SC1 AIC1])
disp(["       St = 2 OLS"])
disp("        Bhat2           Sig2hat2          varbhat          stde                 tval                  pval          R^2          R^2bar           SC            AIC")
disp([bhat2 sig2hat2 varbhat2 stde2 t_val2 p_val2 R22 R2_2 SC2 AIC2])
disp("==============================================================================================")
disp(["       True Parameter"])
disp(["         B1           Sig2_1"])
disp([Beta_0 Sig2_1])
disp(["         B2           Sig2_2"])
disp([Beta_1 Sig2_2])
disp("==============================================================================================")
          