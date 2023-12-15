%% parameter constraints

function [valid] = paramconst(theta,Data)

validm = ones(10,1);

% Sig2 > 0
validm(1) = abs(theta(3)) < 1;   % abs(Phi) < 1
validm(2) = abs(theta(4)) < 1;
validm(3) = abs(theta(5)) < 1;
validm(4) = abs(theta(6)) < 1;

validm(5) = theta(7) > 0;        % Sig^2_1 > 0
validm(6) = theta(7) > 0;        % Sig^2_2 > 0
 
validm(7) = theta(10) < 1;        % 0 < p0 < 1
validm(8) = theta(10) > 0;

validm(9) = theta(9) < 1;        % 0 < q0 < 1
validm(10) = theta(9) > 0;

validm(11) = theta(11) < 1;        % 0 < p1 < 1
validm(12) = theta(11) > 0;

validm(13) = theta(12) < 1;        % 0 < q1 < 1
validm(14) = theta(12) > 0;

valid = min(validm); 

end
