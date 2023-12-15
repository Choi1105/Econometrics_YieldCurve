%% Compute log likelihood function

function [lnL, Spm] = lnlik(theta,Data)

% Data
T = rows(Data);
Y = Data(:,1);
X = Data(:,2:3);

% Parameter
beta1 = theta(1:2);
beta2 = theta(3:4);
beta3 = theta(5:6);
sig21 = theta(7);
sig22 = theta(8);
sig23 = theta(9);

p11 = theta(10);
p22 = theta(11);

ms = 3;
P = zeros(ms, ms);
P(1,1) = p11;
P(1,2) = 1 - p11;
P(1,3) = 0; % jump 불가능.
P(2,1) = 0;
P(2,2) = p22;
P(2,3) = 1 - p22;
P(3,1) = 0;
P(3,2) = 0;
P(3,3) = 1; % 뒤로 불가능.

pr1 = 1;
pr2 = 0;
pr3 = 0;

Spm = zeros(T,3);

lnLm = zeros(T, 1);

for t = 1:T
    prtl1 = pr1*P(1,1) + pr2*P(2,1) + pr3*P(3,1);
    prtl2 = pr1*P(1,2) + pr2*P(2,2) + pr3*P(3,2);
    prtl3 = pr1*P(1,3) + pr2*P(2,3) + pr3*P(3,3);

    pdf1 = normpdf(Y(t), X(t,:)*beta1, sqrt(sig21));
    pdf2 = normpdf(Y(t), X(t,:)*beta2, sqrt(sig22));
    pdf3 = normpdf(Y(t), X(t,:)*beta3, sqrt(sig23));

    pdf = prtl1*pdf1 + prtl2*pdf2 + prtl3*pdf3;

    lnLm(t) = log(pdf);

    prtt1 = prtl1*pdf1/pdf;
    prtt2 = prtl2*pdf2/pdf;
    prtt3 = prtl3*pdf3/pdf;

    Spm(t,:) = [prtt1, prtt2, prtt3];

    pr1 = prtt1;
    pr2 = prtt2;
    pr3 = prtt3;
end

lnL = sum(lnLm);

end
