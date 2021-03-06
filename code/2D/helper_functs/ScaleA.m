function [A,AH,b] = ScaleA(n,A,AH,b)

% Scales mu, A and f so that the largest eigenvalue of A'*A is 1 and the
% new problem
%
% min sum_i (||wi|| + beta/2 ||Diu - wi||^2) + mu/2 ||Au - b||^2
%
% is equivalent to the old one.  
%
% If option is assigned, mu will be scaled accordingly.
%
% Written by: Chengbo Li


% check if A*A'=I
tmp = rand(length(b),1);
if norm(A(AH(tmp))-tmp,1)/norm(tmp,1) < 1e-3
    return;
end

eopts.disp = 0;
eopts.tol = .05;
if ~isreal(A(rand(n,1)))
    eopts.isreal = false;
end

fh = @(x) AH(A(x));
s2 = eigs(fh,n,1,'lm',eopts);
if real(s2) > 1 + 1e-10
    b = b/sqrt(s2);
    A = @(x) A(x)/sqrt(s2);
    AH = @(x) AH(x)/sqrt(s2); 
end

return