function [sigma] = conc_factor(k)
%%%%%%%%%%%%%%%%%%
% Concentration factors
% for this work, exponential factor is used, code originally from
% Theresa Scarnati, 2019


% Previously, the polynomial factor in Gelb 2011 eq(11) was used
% as in paper, p=5 was used
%%%%%%%%%%%%%%%%%%

% polynomial factor
% p=1; % it says 5 above, but 1 works better in code (will use exp anyway later)
sigma = 1i*pi*p*k.^p; %1i*pi*p*k.^p;
sigma = 2 * sigma / sqrt(numel(k)); % fix scaling issues

% exponential factor
% order = 2; % not sure what importance this has...
% eta   = abs(k)/max(abs(k));
% fun   = @(x) exp(1./(order*x.*(x-1)));
% C     = pi/integral(fun,1/max(abs(k),[], 'all'),1-1/max(abs(k),[], 'all'));
% sigma   = C*eta.*exp(1./(order*eta.*(eta-1)));
% sigma(abs(eta-1)<1e-8)=0;


% trig factor
% sigma = sin(pi * k) / sinint(pi);
% sigma = sigma / pi;


end