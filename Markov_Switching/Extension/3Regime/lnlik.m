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
p12 = theta(11);
p21 = theta(12);
p22 = theta(13);
p31 = theta(14);
p32 = theta(15);

ms = 3;

P = zeros(ms);
P(1,1) = p11;
P(1,2) = p12;
P(1,3) = 1 - p11 - p12;
P(2,1) = p21;
P(2,2) = p22;
P(2,3) = 1 - p21 - p22;
P(3,1) = p31;
P(3,2) = p32;
P(3,3) = 1 - p31 - p32;

P_Star = P';

% P_Star(1,1) = p11;
% P_Star(1,2) = p21;
% P_Star(1,3) = p31;
% P_Star(2,1) = p12;
% P_Star(2,2) = p22;
% P_Star(2,3) = p32;
% P_Star(3,1) = P(1,3);
% P_Star(3,2) = P(2,3);
% P_Star(3,3) = P(3,3);

im = [1 1 1]';
IdM = [1 0 0 ; 0 1 0 ; 0 0 1];
A = [IdM - P_Star ; im'];
PROBB = inv(A'*A)*A'*[0;0;0;1];

pr1 = PROBB(1);
pr2 = PROBB(2);
pr3 = PROBB(3);

Spm = zeros(T,ms);

lnLm = zeros(T, 1);

for t = 1:T
    prtl1 = pr1*P(1,1) + pr2*P(2,1) + pr3*P(3,1);
    prtl2 = pr1*P(1,2) + pr2*P(2,2) + pr3*P(3,2);
    prtl3 = pr1*P(1,3) + pr2*P(2,3) + pr3*P(3,3);

    pdf1 = normpdf(Y(t), X(t,:)*beta1, sqrt(sig21));
    pdf2 = normpdf(Y(t), X(t,:)*beta2, sqrt(sig22));
    pdf3 = normpdf(Y(t), X(t,:)*beta3, sqrt(sig23));

    pdf = pdf1*prtl1 + pdf2*prtl2 + pdf3*prtl3;

    prtt1 = pdf1*prtl1/pdf;
    prtt2 = pdf2*prtl2/pdf;
    prtt3 = pdf3*prtl3/pdf;

    lnLm(t) = log(pdf);

    Spm(t,:) = [prtt1, prtt2, prtt3];

    pr1 = prtt1;
    pr2 = prtt2;
    pr3 = prtt3;
end

lnL = sum(lnLm);

end
