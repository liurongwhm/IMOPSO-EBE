function [theta]=SAM(A,B)
AB = A'*B;
Anorm = norm(A,2);
Bnorm = norm(B,2);
V = AB/(Anorm*Bnorm);
theta = acos(V);