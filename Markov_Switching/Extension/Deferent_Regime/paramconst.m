%% parameter constraints

function [valid] = paramconst(theta,Data)

validm = ones(20,1);

% Sig2 > 0
validm(1) = theta(3) > 0;
validm(2) = theta(4) > 0;

validm(3) = theta(5) > 0;
validm(4) = theta(5) < 1;

validm(5) = theta(6) > 0;
validm(6) = theta(6) < 1;

validm(7) = theta(7) > 0;
validm(8) = theta(7) < 1;

validm(9) = theta(8) > 0;
validm(10) = theta(8) < 1;


valid = min(validm); 

end
