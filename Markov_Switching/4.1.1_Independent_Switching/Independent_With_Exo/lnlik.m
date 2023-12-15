%% Compute log likelihood function

function [lnL] = lnlik(theta,Data)

% Data
T = rows(Data);
Y = Data(:,1);
X = Data(:,2);
St = Data(:,3);

% Parameter
beta = theta(1:2); 
sig21 = theta(3);
sig22 = theta(4);
p0 = theta(5);
p1 = theta(6);


% Compute log density for t=1~T
lnLm = zeros(T,1);
for t = 2:T
    St_prob = exp(p0+X(t-1)*p1)/(1+exp(p0+X(t-1)*p1)) * St(t) + (1-exp(p0+X(t-1)*p1)/(1+exp(p0+X(t-1)*p1))) * (1-St(t));
    cond_density = mvnpdf(Y(t),X(t,:)*beta(1)*(1-St(t))+X(t,:)*beta(2)*St(t),sig21*(1-St(t))+sig22*St(t));
    lnLm(t) = log(cond_density * St_prob);

end

% Sum
lnL = sumc(lnLm(2:end));

end