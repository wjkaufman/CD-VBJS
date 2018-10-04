function [j_star,meas_mat] = get_VWJSdata(meas)
% for finding the correct data for final VWJS reconstruction

[size_meas,num_meas] = size(meas); 
l2_comp = zeros(num_meas,num_meas); 
for ii = 1:num_meas 
    for jj = 1:num_meas
        l2_comp(ii,jj) = norm(meas(:,ii)-meas(:,jj),2);
    end
end
tmp = ones(num_meas,1)*max(l2_comp(:)); 
tmp = diag(tmp); 
meas_mat = l2_comp; 
l2_comp = l2_comp + tmp; 

[row,col] = find(l2_comp == min(min(l2_comp))); 
j_star = row(1); 

end

