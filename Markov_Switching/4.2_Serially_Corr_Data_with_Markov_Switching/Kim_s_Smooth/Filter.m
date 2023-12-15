







% Backward (Smoothing equation)

function [Lspm] = Filter(l_p, l_q, Spm, prtlm)
T = rows(Spm);

P = zeros(2,2); % 최종 확률 matrix 4.40 첫번째 분자.
P(1,1) = l_p; 
P(1,2) = 1 - l_p;
P(2,1) = 1 - l_q;
P(2,2) = l_q;

Lspm = zeros(T,2);
Lspm(T,:) = Spm(T,:);

% Backward (Smoothing equation)
t = T-1;
while t >= 1

    % Spm 은 4.37의 결과 3번째 eq 
    prT11 = (P(1,1)*Spm(t,1)*Lspm(t+1,1))/prtlm(t+1,1);  
    prT12 = (P(1,2)*Spm(t,1)*Lspm(t+1,2))/prtlm(t+1,2);
    prT21 = (P(2,1)*Spm(t,2)*Lspm(t+1,1))/prtlm(t+1,1);
    prT22 = (P(2,2)*Spm(t,2)*Lspm(t+1,2))/prtlm(t+1,2);

    prT1 = prT11 + prT12;
    prT2 = prT21 + prT22;

    Lspm(t,:) = [prT1, prT2];

    t = t - 1;

end