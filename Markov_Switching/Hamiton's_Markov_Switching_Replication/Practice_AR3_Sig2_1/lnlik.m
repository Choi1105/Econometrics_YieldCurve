%% Compute log likelihood function

function [lnL, Spm] = lnlik(theta,Data)

% Data
T = rows(Data);
Y = Data(:,1);
X = Data(:,2:4);

% Parameter
b11 = theta(1);
b12 = theta(2);

phi1 = theta(3);
phi2 = theta(4);
phi3 = theta(5);

sig2_1 = theta(6);
sig2_2 = theta(7);

p0 = theta(8);
q0 = theta(9);

P = zeros(2, 2);
P(1,1) = p0;
P(1,2) = 1 - p0;
P(2,1) = 1 - q0;
P(2,2) = q0;

Im = ones(2,1);
Zm = zeros(2,1);

A = [eye(2) - P' ; Im'];
SS_prob = inv(A'*A)*A'*[Zm;1];

% Compute log density for t=1:T
pr1 = SS_prob(1);
pr2 = SS_prob(2);

Spm = zeros(T,2);
prtlm = zeros(T,2);
lnLm = zeros(T, 1);

for t = 1:T
    
    % (4.35)
    prtl1111 = pr1*P(1,1)*P(1,1)*P(1,1);
    prtl1121 = pr1*P(1,2)*P(2,1)*P(1,1);
    prtl1211 = pr1*P(1,1)*P(1,2)*P(2,1);
    prtl2111 = pr1*P(1,1)*P(1,1)*P(1,2);
    prtl2211 = pr1*P(1,1)*P(1,2)*P(2,2);
    prtl1221 = pr1*P(1,2)*P(2,2)*P(2,1);
    prtl2121 = pr1*P(1,2)*P(2,1)*P(1,2);
    prtl2221 = pr1*P(1,2)*P(2,2)*P(2,2);

    prtl2222 = pr2*P(2,2)*P(2,2)*P(2,2);
    prtl2212 = pr2*P(2,1)*P(1,2)*P(2,2);
    prtl2122 = pr2*P(2,2)*P(2,1)*P(1,2);
    prtl1222 = pr2*P(2,2)*P(2,2)*P(2,1);
    prtl1122 = pr2*P(2,2)*P(2,1)*P(1,1);
    prtl1212 = pr2*P(2,1)*P(1,2)*P(2,1);
    prtl2112 = pr2*P(2,1)*P(1,1)*P(1,2);
    prtl1112 = pr2*P(2,1)*P(1,1)*P(1,1);

    prtl1 = prtl1111 + prtl1121 + prtl1211 + prtl1221 + prtl1222 + prtl1122 + prtl1212 + prtl1112;
    prtl2 = prtl2111 + prtl2211 + prtl2121 + prtl2221 + prtl2222 + prtl2212 + prtl2122 + prtl2112;

    prtlm(t,:) = [prtl1, prtl2];

    % (4.33)
    pdf1111 = normpdf(Y(t), b11 + phi1*(X(t,1) - b11) + phi2*(X(t,2)-b11) + phi3*(X(t,3)-b11), sqrt(sig2_1));
    pdf1121 = normpdf(Y(t), b11 + phi1*(X(t,1) - b11) + phi2*(X(t,2)-b12) + phi3*(X(t,3)-b11), sqrt(sig2_1));
    pdf1211 = normpdf(Y(t), b11 + phi1*(X(t,1) - b12) + phi2*(X(t,2)-b11) + phi3*(X(t,3)-b11), sqrt(sig2_1));
    pdf2111 = normpdf(Y(t), b12 + phi1*(X(t,1) - b11) + phi2*(X(t,2)-b11) + phi3*(X(t,3)-b11), sqrt(sig2_2));
    pdf2211 = normpdf(Y(t), b12 + phi1*(X(t,1) - b12) + phi2*(X(t,2)-b11) + phi3*(X(t,3)-b11), sqrt(sig2_2));
    pdf1221 = normpdf(Y(t), b11 + phi1*(X(t,1) - b12) + phi2*(X(t,2)-b12) + phi3*(X(t,3)-b11), sqrt(sig2_1));
    pdf2121 = normpdf(Y(t), b12 + phi1*(X(t,1) - b11) + phi2*(X(t,2)-b12) + phi3*(X(t,3)-b11), sqrt(sig2_2));
    pdf2221 = normpdf(Y(t), b12 + phi1*(X(t,1) - b12) + phi2*(X(t,2)-b12) + phi3*(X(t,3)-b11), sqrt(sig2_2));
    
    pdf2222 = normpdf(Y(t), b12 + phi1*(X(t,1) - b12) + phi2*(X(t,2)-b12) + phi3*(X(t,3)-b12), sqrt(sig2_2));
    pdf2212 = normpdf(Y(t), b12 + phi1*(X(t,1) - b12) + phi2*(X(t,2)-b11) + phi3*(X(t,3)-b12), sqrt(sig2_2));
    pdf2122 = normpdf(Y(t), b12 + phi1*(X(t,1) - b11) + phi2*(X(t,2)-b12) + phi3*(X(t,3)-b12), sqrt(sig2_2));
    pdf1222 = normpdf(Y(t), b11 + phi1*(X(t,1) - b12) + phi2*(X(t,2)-b12) + phi3*(X(t,3)-b12), sqrt(sig2_1));
    pdf1122 = normpdf(Y(t), b11 + phi1*(X(t,1) - b11) + phi2*(X(t,2)-b12) + phi3*(X(t,3)-b12), sqrt(sig2_1));
    pdf1212 = normpdf(Y(t), b11 + phi1*(X(t,1) - b12) + phi2*(X(t,2)-b11) + phi3*(X(t,3)-b12), sqrt(sig2_1));
    pdf2112 = normpdf(Y(t), b12 + phi1*(X(t,1) - b11) + phi2*(X(t,2)-b11) + phi3*(X(t,3)-b12), sqrt(sig2_2));
    pdf1112 = normpdf(Y(t), b11 + phi1*(X(t,1) - b11) + phi2*(X(t,2)-b11) + phi3*(X(t,3)-b12), sqrt(sig2_1));


    pdf = prtl1111*pdf1111 + prtl1121*pdf1121 + prtl1211*pdf1211 + prtl1221*pdf1221 + prtl1222*pdf1222 + prtl1122*pdf1122 + prtl1212*pdf1212 + prtl1112*pdf1112 + ...
          prtl2111*pdf2111 + prtl2211*pdf2211 + prtl2121*pdf2121 + prtl2221*pdf2221 + prtl2222*pdf2222 + prtl2212*pdf2212 + prtl2122*pdf2122 + prtl2112*pdf2112;

    % (4.34)
    lnLm(t) = log(pdf);

    % (4.36)
    prtt1 = (prtl1111*pdf1111 + prtl1121*pdf1121 + prtl1211*pdf1211 + prtl1221*pdf1221 + prtl1222*pdf1222 + prtl1122*pdf1122 + prtl1212*pdf1212 + prtl1112*pdf1112)/pdf;
    prtt2 = (prtl2111*pdf2111 + prtl2211*pdf2211 + prtl2121*pdf2121 + prtl2221*pdf2221 + prtl2222*pdf2222 + prtl2212*pdf2212 + prtl2122*pdf2122 + prtl2112*pdf2112)/pdf;

    Spm(t,:) = [prtt1, prtt2];

    pr1 = prtt1;
    pr2 = prtt2;


lnL = sum(lnLm);

end
