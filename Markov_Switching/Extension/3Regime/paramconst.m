%% parameter constraints

function [valid] = paramconst(theta,Data)

validm = ones(30,1);

% Sig2 > 0
validm(1) = theta(7) > 0;
validm(2) = theta(8) > 0;
validm(3) = theta(9) > 0;

validm(4) = theta(10) > 0;
validm(5) = theta(10) < 1;

validm(6) = theta(11) > 0;
validm(7) = theta(11) < 1;

validm(8) = theta(12) > 0;
validm(9) = theta(12) < 1;

validm(10) = theta(13) > 0;
validm(11) = theta(13) < 1;

validm(12) = theta(14) > 0;
validm(13) = theta(14) < 1;

validm(14) = theta(15) > 0;
validm(15) = theta(15) < 1;

validm(16) = theta(10) + theta(11)  < 1; % 두 확률의 합이 1보다 작아야 한다는 제약
validm(17) = theta(12) + theta(13)  < 1;
validm(18) = theta(14) + theta(15)  < 1;

valid = min(validm); 

end
