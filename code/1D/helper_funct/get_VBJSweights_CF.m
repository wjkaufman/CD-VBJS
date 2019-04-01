function [w,v] = get_VBJSweights_CF(f_tild)
% get VWJS weights
%% CHANGE IN SLOPE NOT CHANGE IN SIGN! 

[N,~] = size(f_tild); 
v = var(f_tild,[],2); 
S = minmod(f_tild);

if sum(S) == 0 
    w = 1/(v+1e-6); 
else
    
    T = abs(S.*v)./max(abs(S.*v));
    ind = find(T<=1/N);
    c = length(find(T>1/N));
    
    % w = zeros(N,1);
    w = 1 - T;
    w(ind) = c.*w(ind);

end


% tmp = 1; 

end

