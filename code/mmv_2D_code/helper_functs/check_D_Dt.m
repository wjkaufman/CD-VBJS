function [flg,x,y] = check_D_Dt(D,Dt,N)

% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 08/24/2016


% Check if Dt is the adjoint of D
% If it is the true adjoint, then we should have x=y
% flg returns true if Dt is true adjoint and false otherwise

u = rand(N) + 1i*rand(N);
theta = 1; 
Du = D(u,theta);
% Du = Du(:,2); 

v = rand(size(Du)) + 1i*rand(size(Du));

Dtv = Dt(v,theta);

x = sum(Du(:).*conj(v(:)));
y = sum(u(:).*conj(Dtv(:)));

tol = 1e-3;

flg = abs(abs(x)-abs(y))/abs(x)<tol;

