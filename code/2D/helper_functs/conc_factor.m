function [sigma] = conc_factor(k)
%%%%%%%%%%%%%%%%%%
% Concentration factors
% for this work, the polynomial factor in Gelb 2011 eq(11) is used
% as in paper, p=5 is used
%%%%%%%%%%%%%%%%%%
p=1; % it says 5 above, but 1 works better in code (will use exp anyway later)
sigma = 1i*pi*p*k.^p; %1i*pi*p*k.^p;
sigma = 2 * sigma / sqrt(numel(k)); % fix scaling issues
end