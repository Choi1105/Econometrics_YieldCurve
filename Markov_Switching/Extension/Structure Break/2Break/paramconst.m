%% parameter constraints

function [valid] = paramconst(theta,Data)

validm = ones(10,1);

% Sig2 > 0
validm(2) = theta(5) > 0;
validm(3) = theta(6) > 0;
validm(4) = theta(7) < 1;
validm(5) = theta(7) > 0;

valid = min(validm); 

end
