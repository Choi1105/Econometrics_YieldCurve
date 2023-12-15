%% parameter constraints

function [valid] = paramconst(theta,Data)

validm = ones(10,1);

% Sig2 > 0
validm(1) = theta(7) > 0;
validm(2) = theta(8) > 0;
validm(3) = theta(9) > 0;
validm(4) = theta(10) > 0;
validm(5) = theta(10) < 1;
validm(6) = theta(11) > 0;
validm(7) = theta(11) < 1;

valid = min(validm); 

end
