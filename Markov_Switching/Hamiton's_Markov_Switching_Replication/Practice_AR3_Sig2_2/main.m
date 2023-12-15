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
DATA = readmatrix("Data.xlsx" , "Range", "C3:C254");

figure
parcorr(DATA);
hY3 = adftest(DATA,Model="ts",Lags=4);

Ym_New = DATA;
nT = rows(Ym_New);

%% Step 2: Estimation %%
% Data
Y = Ym_New(4:end); 
X = [Ym_New(1:nT-3) Ym_New(2:nT-2) Ym_New(3:nT-1)]; 
Data = [Y X];
T = rows(X);
k = cols(X);

% initial value 
% Beta1 Beta2 Phi1 Phi2 Phi3 Sig2_1 Sig2_2 p0 q0
theta0 = [0.4 ; 2.6 ; 0.16 ;0.1; 0.2 ; 0.3 ; 0.1 ; 0.8 ; 0.8];

% index
index = [1;2;3;4;5;6;7;8];
printi = 1;

% Optimization
[thetamx, fmax, Vj, Vinv] = SA_Newton(@lnlik,@paramconst,theta0,Data,printi,index);
[~, Spm] = lnlik(thetamx, Data);


%% Step 3: Results %%
Data = DATA(4:end,1);

% Plot
figure;
yyaxis left;  % 왼쪽 y축 사용
plot(Data, 'color', [0, 0, 1], 'LineWidth', 1);
ylabel('Data', 'FontSize', 12);
ylim([min(Data)-2 max(Data)+2]);  % 데이터에 맞게 y축 범위 설정

yyaxis right;  % 오른쪽 y축 사용
plot(Spm(:, 1), 'r', 'LineWidth', 1);
ylabel('Spm(:,1)', 'FontSize', 12);
ylim([min(Spm(:, 1)) max(Spm(:, 1))+0.1]);  % 데이터에 맞게 y축 범위 설정

% 레전드 추가
legend('Data', 'Spm(:,1)', 'FontSize', 12);