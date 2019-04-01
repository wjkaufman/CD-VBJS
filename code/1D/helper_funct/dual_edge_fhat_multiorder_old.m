%% Dual space edge detector
function [f_jump,sig] = dual_edge_fhat_multiorder_old(order,N,f_hat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f_jump = DualEdge(factor_type,N,f)
% factor_type: (1) trigonometric (2) polynomial (3) exponential
%           N: number of grid points (N) fourier coeff
%           f: function to find the edge map (2*pi periodic)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k     = [0:N/2-1 N/2 -N/2+1:-1]';

% concentration factor input
eta   = abs(k)/max(abs(k));
fun   = @(x) exp(1./(order*x.*(x-1)));
C     = pi/quad(fun,1/max(abs(k)),1-1/max(abs(k)));
sig   = C*eta.*exp(1./(order*eta.*(eta-1)));
sig(abs(eta-1)<1e-8)=0;

% jump_kern = exp(1i*x*k');
% coeff_jump = f_hat.*1i.*sign(k).*sig;
% S_N = jump_kern*coeff_jump; 
% f_jump = real(S_N);


C_N = (1i*sign(k).*sig);
S_N = ifft(C_N.*f_hat); % dual space edge detector
f_jump = real(S_N);