%% Compute log likelihood function

function [lnL, Spm] = lnlik(theta,Data)

% Data
T = rows(Data);
Y = Data(:,1);
X = Data(:,2:3);

% Parameter
b11 = theta(1);
b12 = theta(2);
phi1 = theta(3);
sig21 = theta(4);
sig22 = theta(5);
p0 = theta(6);
q0 = theta(7);

pr1 = (1 - p0) / (2 - p0 - q0);
pr2 = (1 - q0) / (2 - p0 - q0);

P_Star = zeros(2,2);
P_Star(1,1) = p0;
P_Star(1,2) = 1 - p0;
P_Star(2,1) = 1 - q0;
P_Star(2,2) = q0;

im = [1 1]';

PiMat = [pr1;pr2];

IdM = [1 0 ; 0 1];

(IdM - P_Star)*PiMat

A = [IdM - P_Star ; im'];

inv(A'*A)*A'*[0;0;1];


Spm = zeros(T,2);

lnLm = zeros(T, 1);

for t = 1:T
    prtl11 = pr1*P_Star(1,1);
    prtl12 = pr1*P_Star(1,2);
    prtl21 = pr2*P_Star(2,1);
    prtl22 = pr2*P_Star(2,2);

    pdf11 = normpdf(Y(t), b11 + phi1*(X(t) - b11), sqrt(sig21));
    pdf12 = normpdf(Y(t), b12 + phi1*(X(t) - b11), sqrt(sig22));
    pdf21 = normpdf(Y(t), b11 + phi1*(X(t) - b12), sqrt(sig21));
    pdf22 = normpdf(Y(t), b12 + phi1*(X(t) - b12), sqrt(sig22));
    
    pdf = prtl11*pdf11 + prtl12*pdf12 + prtl21*pdf21 + prtl22*pdf22;
    
    lnLm(t) = log(pdf);
    
    prtt11 = prtl11*pdf11/pdf;
    prtt12 = prtl12*pdf12/pdf;
    prtt21 = prtl21*pdf21/pdf;
    prtt22 = prtl22*pdf22/pdf;
    prtI1 = prtt11 + prtt21;
    prtI2 = prtt12 + prtt22;
    Spm(t,:) = [prtI1, prtI2];
    pr1 = prtI1;
    pr2 = prtI2;
end

lnL = sum(lnLm);

end
