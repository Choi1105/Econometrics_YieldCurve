p0 = 0.95;
q0 = 0.95;

pr1 = (1 - p0) / (2 - p0 - q0);
pr2 = (1 - q0) / (2 - p0 - q0);

P_Star = zeros(2,2);
P_Star(1,1) = p0;
P_Star(1,2) = 1 - p0;
P_Star(2,1) = 1 - q0;
P_Star(2,2) = q0;

im = [1 1]';

PiMat = [pr1;pr2];

IdM = [1 0 ; 0 1];

(IdM - P_Star)*PiMat

A = [IdM - P_Star ; im'];

A*PiMat

inv(A'*A)*A'*[0;0;1]
