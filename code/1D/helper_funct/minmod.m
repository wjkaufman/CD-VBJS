function [mm_f] = minmod(f)
% here, each column of f is a measurement, that is f is N x J with J
% measurements each having length N. 

[N,J] = size(f); 
signs = sign(f); 
diff_signs = diff(signs,1,2);
sums = sum(abs(diff_signs),2); 
ind = find(sums == 0); 
mm_f = zeros(N,1); 
mm_f(ind) = signs(ind,1).*min(abs(f(ind,:)),[],2);
    

end

