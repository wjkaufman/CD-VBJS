function [mm_f] = minmod(f)
% f is assumed to be a multi-dimensional array, where
% the first D-1 dimensions correspond to a single measurement, and
% the last dimension corresponds to different realizations of f
% here, each column of f is a measurement, that is f is N x J with J
% measurements each having length N. 

f_size = size(f); 
signs = sign(f);
diff_signs = diff(signs,1,length(f_size));
sums = sum(abs(diff_signs),length(f_size));
mm_f = (sums == 0) .* mean(signs,length(f_size)).*min(abs((sums == 0).*f), ...
            [],length(f_size));

end


% (ind + ...
            % includes all 
            %(0:(f_size(end)-1))*prod(f_size(1:(end-1))))