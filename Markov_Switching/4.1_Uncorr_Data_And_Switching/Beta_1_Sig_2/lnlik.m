%% Compute log likelihood function

function [lnL] = lnlik(theta,Data)

% Data
T = rows(Data);
Y = Data(:,1);
X = Data(:,2:end);

St = zeros(T,1);
for t = [1 5 8]
    St(t*1000:(t*1000)+1000,1) = 1;
end

% Parameter
beta = theta(1); 
sig21 = theta(2);
sig22 = theta(3);

% Compute log density for t=1~T
lnLm = zeros(T,1);
for t = 1:T

    lnLm(t) = log(mvnpdf(Y(t),X(t,:)*beta(1)*(1-St(t))+X(t,:)*beta(1)*St(t),sig21*(1-St(t))+sig22*St(t)));

end

% Sum
lnL = sumc(lnLm);

end