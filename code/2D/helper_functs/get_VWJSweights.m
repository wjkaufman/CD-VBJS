function [w,v] = get_VWJSweights(sparse_meas,tau)
% get VWJS weights

v = var(sparse_meas,1,2); 
norm_F = abs(sparse_meas)./max(abs(sparse_meas),[],1); 
% figure; plot(norm_F); 
% hold on ; 
% plot([1,length(norm_F)],[tau,tau],'k--'); 

I = find(mean(norm_F,2)>tau); 
C = mean(sum(abs(sparse_meas),1)); 
while C < 1 
    C = 10*C;
end
w = C*(1-v/max(v));
w(I) = 1/C*(1-v(I)/max(v(I)));

% if C < 1
%     w = 1/C*(1-v/max(v));
%     w(I) = C*(1-v(I)/max(v(I)));
% else
%     w = C*(1-v/max(v));
%     w(I) = 1/C*(1-v(I)/max(v(I)));
% end

end

