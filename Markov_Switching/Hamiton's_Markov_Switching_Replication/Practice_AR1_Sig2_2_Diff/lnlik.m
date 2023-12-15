%% Compute log likelihood function

function [lnL, Spm, Rpm] = lnlik(theta,Data)

% Data
T = rows(Data);
Y = Data(:,1);
X = Data(:,2:5);

% Parameter
b11 = theta(1);
b12 = theta(2);

phi1 = theta(3);
phi2 = theta(4);
phi3 = theta(5);
phi4 = theta(6);

sig2_1 = theta(7);
sig2_2 = theta(8);

p0 = theta(9);
q0 = theta(10);

p1 = theta(11);
q1 = theta(12);

P = zeros(2,2);
P(1,1) = p0;
P(1,2) = 1 - p0;
P(2,1) = 1 - q0;
P(2,2) = q0;

PPAP = zeros(2,2);
PPAP(1,1) = p1;
PPAP(1,2) = 1 - p1;
PPAP(2,1) = 1 - q1;
PPAP(2,2) = q1;

pr11 = (1 - p1) / (2 - p1 - q1);
pr22 = (1 - q1) / (2 - p1 - q1);

Im = ones(2,1);
Zm = zeros(2,1);

A = [eye(2) - P' ; Im'];
SS_prob = inv(A'*A)*A'*[Zm;1];

% Compute log density for t=1:T
pr1 = SS_prob(1);
pr2 = SS_prob(2);

Spm = zeros(T,2);
Rpm = zeros(T,2);
%prtlm = zeros(T,2);
lnLm = zeros(T, 1);

for t = 1:T
    prtl11 = pr1*P(1,1);
    prtl12 = pr1*P(1,2);
    prtl21 = pr2*P(2,1);
    prtl22 = pr2*P(2,2);

    prtlV11 = pr11*PPAP(1,1);
    prtlV12 = pr11*PPAP(1,2);
    prtlV21 = pr22*PPAP(2,1);
    prtlV22 = pr22*PPAP(2,2);

    prtls1 = prtl11 + prtl21;
    prtls2 = prtl12 + prtl22;
    prtlr1 = prtlV11 + prtlV21;
    prtlr2 = prtlV12 + prtlV22;

    pdf11 = normpdf(Y(t), b11 + phi1*X(t,1) + phi2*X(t,2) + phi3*X(t,3) + phi4*X(t,4), sqrt(sig2_1));
    pdf12 = normpdf(Y(t), b12 + phi1*X(t,1) + phi2*X(t,2) + phi3*X(t,3) + phi4*X(t,4), sqrt(sig2_2));
    pdf21 = normpdf(Y(t), b11 + phi1*X(t,1) + phi2*X(t,2) + phi3*X(t,3) + phi4*X(t,4), sqrt(sig2_2));
    pdf22 = normpdf(Y(t), b12 + phi1*X(t,1) + phi2*X(t,2) + phi3*X(t,3) + phi4*X(t,4), sqrt(sig2_1));

    pdf = pdf11*prtls1*prtlr1 + pdf21*prtls2*prtlr1 + pdf12*prtls1*prtlr2 + pdf22*prtls2*prtlr2;
    
    lnLm(t) = log(pdf);
    
    lnLm(t) = log(pdf);
    
    prtts1 = (pdf11*prtls1*prtlr1 + pdf12*prtls1*prtlr2)/pdf;
    prtts2 = (pdf21*prtls2*prtlr1 + pdf22*prtls2*prtlr2)/pdf;
    prttr1 = (pdf11*prtls1*prtlr1 + pdf21*prtls2*prtlr1)/pdf;
    prttr2 = (pdf12*prtls1*prtlr2 + pdf22*prtls2*prtlr2)/pdf;
    
    Spm(t,:) = [prtts1, prtts2];
    Rpm(t,:) = [prttr1, prttr2];


    pr1 = prtts1;
    pr2 = prtts2;
    pr11 = prttr1;
    pr22 = prttr2;
end

lnL = sum(lnLm);

end
