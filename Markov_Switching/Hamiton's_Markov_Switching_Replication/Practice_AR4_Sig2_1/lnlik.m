%% Compute log likelihood function

function [lnL, Spm] = lnlik(theta,Data)

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

sig2 = theta(7);

q0 = theta(8);
p0 = theta(9);

P = zeros(2, 2);
P(1,1) = q0;
P(1,2) = 1 - q0;
P(2,1) = 1 - p0;
P(2,2) = p0;

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
    prtl11111 = pr1*P(1,1)*P(1,1)*P(1,1)*P(1,1);
    prtl11121 = pr1*P(1,2)*P(2,1)*P(1,1)*P(1,1);
    prtl11211 = pr1*P(1,1)*P(1,2)*P(2,1)*P(1,1);
    prtl12111 = pr1*P(1,1)*P(1,1)*P(1,2)*P(2,1);
    prtl21111 = pr1*P(1,1)*P(1,1)*P(1,1)*P(1,2);
    prtl22111 = pr1*P(1,1)*P(1,1)*P(1,2)*P(2,2);
    prtl21211 = pr1*P(1,1)*P(1,2)*P(2,1)*P(1,2);
    prtl21121 = pr1*P(1,2)*P(2,1)*P(1,1)*P(1,2);
    prtl12211 = pr1*P(1,1)*P(1,2)*P(2,2)*P(2,1);
    prtl12121 = pr1*P(1,2)*P(2,1)*P(1,2)*P(2,1);
    prtl11221 = pr1*P(1,2)*P(2,2)*P(2,1)*P(1,1);
    prtl22211 = pr1*P(1,1)*P(1,2)*P(2,2)*P(2,2);
    prtl22121 = pr1*P(1,2)*P(2,1)*P(1,2)*P(2,2);
    prtl21221 = pr1*P(1,2)*P(2,2)*P(2,1)*P(1,2);
    prtl12221 = pr1*P(1,2)*P(2,2)*P(2,2)*P(1,2);
    prtl22221 = pr1*P(1,2)*P(2,2)*P(2,2)*P(2,2);

    prtl22222 = pr2*P(2,2)*P(2,2)*P(2,2)*P(2,2);
    prtl22212 = pr2*P(2,1)*P(1,2)*P(2,2)*P(2,2);
    prtl22122 = pr2*P(2,2)*P(2,1)*P(1,2)*P(2,2);
    prtl21222 = pr2*P(2,2)*P(2,2)*P(2,1)*P(1,2);
    prtl12222 = pr2*P(2,2)*P(2,2)*P(2,2)*P(2,1);
    prtl11222 = pr2*P(2,2)*P(2,2)*P(2,1)*P(1,1);
    prtl12122 = pr2*P(2,2)*P(2,1)*P(1,2)*P(2,1);
    prtl12212 = pr2*P(2,1)*P(1,2)*P(2,2)*P(2,1);
    prtl21122 = pr2*P(2,2)*P(2,1)*P(1,1)*P(1,2);
    prtl21212 = pr2*P(2,1)*P(1,2)*P(2,1)*P(1,2);
    prtl22112 = pr2*P(2,1)*P(1,1)*P(1,2)*P(2,2);
    prtl11122 = pr2*P(2,2)*P(2,1)*P(1,1)*P(1,1);
    prtl11212 = pr2*P(2,1)*P(1,2)*P(2,1)*P(1,1);
    prtl12112 = pr2*P(2,1)*P(1,1)*P(1,2)*P(2,1);
    prtl21112 = pr2*P(2,1)*P(1,1)*P(1,1)*P(1,2);
    prtl11112 = pr2*P(2,1)*P(1,1)*P(1,1)*P(1,1);

    prtl1 = prtl11111 + prtl11121 + prtl11211 + prtl12111 + prtl12211 + prtl12121 + prtl11221 + prtl12221 + ...
            prtl12222 + prtl11222 + prtl12122 + prtl12212 + prtl11122 + prtl11212 + prtl12112 + prtl11112;
    prtl2 = prtl21111 + prtl22111 + prtl21211 + prtl21121 + prtl22211 + prtl22121 + prtl21221 + prtl22221 + ...
            prtl22222 + prtl22212 + prtl22122 + prtl21222 + prtl21122 + prtl21212 + prtl22112 + prtl21112;

    prtlm(t,:) = [prtl1, prtl2];

    % (4.33)
    pdf11111 = normpdf(Y(t), b11 + phi1*(X(t,1) - b11) + phi2*(X(t,2)-b11) + phi3*(X(t,3)-b11) + phi4*(X(t,4) - b11), sqrt(sig2));
    pdf11121 = normpdf(Y(t), b11 + phi1*(X(t,1) - b11) + phi2*(X(t,2)-b11) + phi3*(X(t,3)-b12) + phi4*(X(t,4) - b11), sqrt(sig2));
    pdf11211 = normpdf(Y(t), b11 + phi1*(X(t,1) - b11) + phi2*(X(t,2)-b12) + phi3*(X(t,3)-b11) + phi4*(X(t,4) - b11), sqrt(sig2));
    pdf12111 = normpdf(Y(t), b11 + phi1*(X(t,1) - b12) + phi2*(X(t,2)-b11) + phi3*(X(t,3)-b11) + phi4*(X(t,4) - b11), sqrt(sig2));
    pdf12211 = normpdf(Y(t), b11 + phi1*(X(t,1) - b12) + phi2*(X(t,2)-b12) + phi3*(X(t,3)-b11) + phi4*(X(t,4) - b11), sqrt(sig2));
    pdf12121 = normpdf(Y(t), b11 + phi1*(X(t,1) - b12) + phi2*(X(t,2)-b11) + phi3*(X(t,3)-b12) + phi4*(X(t,4) - b11), sqrt(sig2));
    pdf11221 = normpdf(Y(t), b11 + phi1*(X(t,1) - b11) + phi2*(X(t,2)-b12) + phi3*(X(t,3)-b12) + phi4*(X(t,4) - b11), sqrt(sig2));
    pdf12221 = normpdf(Y(t), b11 + phi1*(X(t,1) - b12) + phi2*(X(t,2)-b12) + phi3*(X(t,3)-b12) + phi4*(X(t,4) - b11), sqrt(sig2));
    pdf12222 = normpdf(Y(t), b11 + phi1*(X(t,1) - b12) + phi2*(X(t,2)-b12) + phi3*(X(t,3)-b12) + phi4*(X(t,4) - b12), sqrt(sig2));
    pdf11222 = normpdf(Y(t), b11 + phi1*(X(t,1) - b11) + phi2*(X(t,2)-b12) + phi3*(X(t,3)-b12) + phi4*(X(t,4) - b12), sqrt(sig2));
    pdf12122 = normpdf(Y(t), b11 + phi1*(X(t,1) - b12) + phi2*(X(t,2)-b11) + phi3*(X(t,3)-b12) + phi4*(X(t,4) - b12), sqrt(sig2));
    pdf12212 = normpdf(Y(t), b11 + phi1*(X(t,1) - b11) + phi2*(X(t,2)-b11) + phi3*(X(t,3)-b11) + phi4*(X(t,4) - b11), sqrt(sig2));
    pdf11122 = normpdf(Y(t), b11 + phi1*(X(t,1) - b11) + phi2*(X(t,2)-b11) + phi3*(X(t,3)-b12) + phi4*(X(t,4) - b12), sqrt(sig2));
    pdf11212 = normpdf(Y(t), b11 + phi1*(X(t,1) - b11) + phi2*(X(t,2)-b12) + phi3*(X(t,3)-b11) + phi4*(X(t,4) - b12), sqrt(sig2));
    pdf12112 = normpdf(Y(t), b11 + phi1*(X(t,1) - b12) + phi2*(X(t,2)-b11) + phi3*(X(t,3)-b11) + phi4*(X(t,4) - b12), sqrt(sig2));
    pdf11112 = normpdf(Y(t), b11 + phi1*(X(t,1) - b11) + phi2*(X(t,2)-b11) + phi3*(X(t,3)-b11) + phi4*(X(t,4) - b12), sqrt(sig2));

    pdf21111 = normpdf(Y(t), b12 + phi1*(X(t,1) - b11) + phi2*(X(t,2)-b11) + phi3*(X(t,3)-b11) + phi4*(X(t,4) - b11), sqrt(sig2));
    pdf22111 = normpdf(Y(t), b12 + phi1*(X(t,1) - b12) + phi2*(X(t,2)-b11) + phi3*(X(t,3)-b11) + phi4*(X(t,4) - b11), sqrt(sig2));
    pdf21211 = normpdf(Y(t), b12 + phi1*(X(t,1) - b11) + phi2*(X(t,2)-b12) + phi3*(X(t,3)-b11) + phi4*(X(t,4) - b11), sqrt(sig2));
    pdf21121 = normpdf(Y(t), b12 + phi1*(X(t,1) - b11) + phi2*(X(t,2)-b11) + phi3*(X(t,3)-b12) + phi4*(X(t,4) - b11), sqrt(sig2));
    pdf22211 = normpdf(Y(t), b12 + phi1*(X(t,1) - b12) + phi2*(X(t,2)-b12) + phi3*(X(t,3)-b11) + phi4*(X(t,4) - b11), sqrt(sig2));
    pdf22121 = normpdf(Y(t), b12 + phi1*(X(t,1) - b12) + phi2*(X(t,2)-b11) + phi3*(X(t,3)-b12) + phi4*(X(t,4) - b11), sqrt(sig2));
    pdf21221 = normpdf(Y(t), b12 + phi1*(X(t,1) - b11) + phi2*(X(t,2)-b12) + phi3*(X(t,3)-b12) + phi4*(X(t,4) - b11), sqrt(sig2));
    pdf22221 = normpdf(Y(t), b12 + phi1*(X(t,1) - b12) + phi2*(X(t,2)-b12) + phi3*(X(t,3)-b12) + phi4*(X(t,4) - b11), sqrt(sig2));
    pdf22222 = normpdf(Y(t), b12 + phi1*(X(t,1) - b12) + phi2*(X(t,2)-b12) + phi3*(X(t,3)-b12) + phi4*(X(t,4) - b12), sqrt(sig2));
    pdf22212 = normpdf(Y(t), b12 + phi1*(X(t,1) - b12) + phi2*(X(t,2)-b12) + phi3*(X(t,3)-b11) + phi4*(X(t,4) - b12), sqrt(sig2));
    pdf22122 = normpdf(Y(t), b12 + phi1*(X(t,1) - b12) + phi2*(X(t,2)-b11) + phi3*(X(t,3)-b12) + phi4*(X(t,4) - b12), sqrt(sig2));
    pdf21222 = normpdf(Y(t), b12 + phi1*(X(t,1) - b11) + phi2*(X(t,2)-b12) + phi3*(X(t,3)-b12) + phi4*(X(t,4) - b12), sqrt(sig2));
    pdf21122 = normpdf(Y(t), b12 + phi1*(X(t,1) - b11) + phi2*(X(t,2)-b11) + phi3*(X(t,3)-b12) + phi4*(X(t,4) - b12), sqrt(sig2));
    pdf21212 = normpdf(Y(t), b12 + phi1*(X(t,1) - b11) + phi2*(X(t,2)-b12) + phi3*(X(t,3)-b11) + phi4*(X(t,4) - b12), sqrt(sig2));
    pdf22112 = normpdf(Y(t), b12 + phi1*(X(t,1) - b12) + phi2*(X(t,2)-b11) + phi3*(X(t,3)-b11) + phi4*(X(t,4) - b12), sqrt(sig2));
    pdf21112 = normpdf(Y(t), b12 + phi1*(X(t,1) - b11) + phi2*(X(t,2)-b11) + phi3*(X(t,3)-b11) + phi4*(X(t,4) - b12), sqrt(sig2));

    pdf = prtl11111*pdf11111 + prtl11121*pdf11121 + prtl11211*pdf11211 + prtl12111*pdf12111 + prtl12211*pdf12211 + prtl12121*pdf12121 + prtl11221*pdf11221 + prtl12221*pdf12221 + ...
          prtl12222*pdf12222 + prtl11222*pdf11222 + prtl12122*pdf12122 + prtl12212*pdf12212 + prtl11122*pdf11122 + prtl11212*pdf11212 + prtl12112*pdf12112 + prtl11112*pdf11112 + ...
          prtl21111*pdf21111 + prtl22111*pdf22111 + prtl21211*pdf21211 + prtl21121*pdf21121 + prtl22211*pdf22211 + prtl22121*pdf22121 + prtl21221*pdf21221 + prtl22221*pdf22221 + ...
          prtl22222*pdf22222 + prtl22212*pdf22212 + prtl22122*pdf22122 + prtl21222*pdf21222 + prtl21122*pdf21122 + prtl21212*pdf21212 + prtl22112*pdf22112 + prtl21112*pdf21112;

    % (4.34)
    lnLm(t) = log(pdf);

    % (4.36)
    prtt1 = (prtl11111*pdf11111 + prtl11121*pdf11121 + prtl11211*pdf11211 + prtl12111*pdf12111 + prtl12211*pdf12211 + prtl12121*pdf12121 + prtl11221*pdf11221 + prtl12221*pdf12221 + ...
          prtl12222*pdf12222 + prtl11222*pdf11222 + prtl12122*pdf12122 + prtl12212*pdf12212 + prtl11122*pdf11122 + prtl11212*pdf11212 + prtl12112*pdf12112 + prtl11112*pdf11112)/pdf;
    prtt2 = (prtl21111*pdf21111 + prtl22111*pdf22111 + prtl21211*pdf21211 + prtl21121*pdf21121 + prtl22211*pdf22211 + prtl22121*pdf22121 + prtl21221*pdf21221 + prtl22221*pdf22221 + ...
          prtl22222*pdf22222 + prtl22212*pdf22212 + prtl22122*pdf22122 + prtl21222*pdf21222 + prtl21122*pdf21122 + prtl21212*pdf21212 + prtl22112*pdf22112 + prtl21112*pdf21112)/pdf;

    Spm(t,:) = [prtt1, prtt2];

    pr1 = prtt1;
    pr2 = prtt2;



lnL = sum(lnLm);

end
