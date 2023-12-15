%% Compute log likelihood function

function [lnL] = lnlik(theta,Data)

% Data
T = rows(Data);
Y = Data(:,1);
X = Data(:,2);

% Parameter
beta = theta(1:2); 
sig21 = theta(3);
sig22 = theta(4);
p0 = theta(5);

% Compute log density for t=1~T
lnLm = zeros(T,1);
for t = 1:T
    cond_density1 = normpdf(Y(t),X(t,:)*beta(1),sig21);
    cond_density2 = normpdf(Y(t),X(t,:)*beta(2),sig22);
    St_prob1 = exp(p0)/(1+exp(p0));
    St_prob2 = 1-((exp(p0))/(1+exp(p0)));
    lnLm(t) = log(cond_density1 * St_prob1 + cond_density2 * St_prob2);
end

% Sum
lnL = sumc(lnLm);

end