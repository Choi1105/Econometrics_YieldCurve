%% Compute log likelihood function

function [lnL, Spm] = lnlik(theta,Data)

% Data
T = rows(Data);
Y = Data(:,1);
X = Data(:,2:3);

% Parameter
beta1 = theta(1:2);
beta2 = theta(3:4);

sig21 = theta(5);
sig22 = theta(6);

p11 = theta(7);

P = zeros(2,2);
P(1,1) = p11;
P(1,2) = 1 - p11;
P(2,1) = 0;
P(2,2) = 1;

pr1 = 1;
pr2 = 0;

Spm = zeros(T,2);

lnLm = zeros(T, 1);

for t = 1:T
    prtl1 = pr1*P(1,1) + pr2*P(2,1);
    prtl2 = pr1*P(1,2) + pr2*P(2,2);

    pdf1 = normpdf(Y(t), X(t,:)*beta1, sqrt(sig21));
    pdf2 = normpdf(Y(t), X(t,:)*beta2, sqrt(sig22));

    pdf = prtl1*pdf1 + prtl2*pdf2;

    lnLm(t) = log(pdf);

    prtt1 = prtl1*pdf1/pdf;
    prtt2 = prtl2*pdf2/pdf;

    Spm(t,:) = [prtt1, prtt2];

    pr1 = prtt1;
    pr2 = prtt2;
end

lnL = sum(lnLm);

end
