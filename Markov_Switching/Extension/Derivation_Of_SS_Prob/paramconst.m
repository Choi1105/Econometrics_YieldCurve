%% parameter constraints

function [valid] = paramconst(theta,Data)

validm = ones(10,1);

% Sig2 > 0
validm(1) = theta(4) > 0;
validm(2) = theta(5) > 0;
validm(3) = theta(6) > 0;
validm(4) = theta(6) < 1;
validm(5) = theta(7) > 0;
validm(6) = theta(7) < 1;
validm(7) = theta(3) < 1;
validm(8) = theta(3) > -1;


valid = min(validm); 

end
