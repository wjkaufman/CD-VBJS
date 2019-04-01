function [j_star,meas_mat] = get_VBJSdata_emd(f_tild)

[N,J] = size(f_tild); 
emd_comp = zeros(J,J); 
nbins = 10; 

for ii = 1:J
    for jj = 1:J
        
        
        [w1, f1] = hist(f_tild(:,ii),nbins); 
        [w2, f2] = hist(f_tild(:,jj),nbins); 
        
        w1 = w1/sum(w1); 
        w2 = w2/sum(w2); 
        
        [~,emd_comp(ii,jj)] = emd(f1',f2',w1',w2', @gdf);
    end 
end
        
tmp = ones(J,1)*max(emd_comp(:)); 
tmp = diag(tmp); 
meas_mat = emd_comp; 
emd_comp = emd_comp + tmp; 

[row,col] = find(emd_comp == min(min(emd_comp))); 
j_star = row(1); 


end

