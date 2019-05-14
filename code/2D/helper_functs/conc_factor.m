function [sigma] = conc_factor(k, order)
%%%%%%%%%%%%%%%%%%
% Concentration factors
% for this work, exponential factor is used, code originally from
% Theresa Scarnati, 2019
%%%%%%%%%%%%%%%%%%

% polynomial factor
% p=1; % it says 5 above, but 1 works better in code (will use exp anyway later)
% sigma = 1i*pi*p*k.^p; %1i*pi*p*k.^p;
% sigma = 2 * sigma / sqrt(numel(k)); % fix scaling issues

% exponential factor
% order: order of the exponential factor (controls spread of sigma)
maxk = max(abs(k), [], 'all');
eta   = abs(k) ./ maxk; % changed from abs(k) ./ ...
fun   = @(x) exp(1./(order*x.*(x-1)));
C     = pi/integral(fun,1/maxk,1-1/maxk);
sigma   = C*eta.*exp(1./(order*eta.*(eta-1)));
sigma(abs(eta-1)<1e-8)=0;
sigma = 1i .* sign(k) .* sigma; % my own twist, trying to replicate good results from poly factor

% trig factor
% sigma = sin(pi * k) / sinint(pi);
% sigma = sigma / pi;
end