function [j_star,meas_mat] = get_VBJSdata_mh(f_tild)

[N,J] = size(f_tild); 
mh_comp = zeros(J,J); 

for ii = 1:J
    for jj = 1:J
        
        [mh_comp(ii,jj)] = ModHausdorffDist(f_tild(:,ii),f_tild(:,jj));
        
    end 
end
        
tmp = ones(J,1)*max(mh_comp(:)); 
tmp = diag(tmp); 
meas_mat = mh_comp; 
mh_comp = mh_comp + tmp; 

[row,col] = find(mh_comp == min(min(mh_comp))); 
j_star = row(1); 


end

