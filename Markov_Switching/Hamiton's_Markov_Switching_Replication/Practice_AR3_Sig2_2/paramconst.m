%% parameter constraints

function [valid] = paramconst(theta,Data)

validm = ones(10,1);

% Sig2 > 0
validm(1) = abs(theta(3)) < 1;   % abs(Phi) < 1
validm(2) = abs(theta(4)) < 1;
validm(3) = abs(theta(5)) < 1;

validm(4) = theta(6) > 0;        % Sig^2 > 0

validm(5) = theta(7) < 1;        % 0 < p0 < 1
validm(6) = theta(7) > 0;

validm(7) = theta(8) < 1;        % 0 < q0 < 1
validm(8) = theta(8) > 0;

valid = min(validm); 

end
