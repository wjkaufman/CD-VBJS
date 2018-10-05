function [D,Dt] = get_D_Dt2(k,p,q,r,opts)

if opts.nonuniform_grid
    [D,Dt] = FD3D_NU2(k,p,q,r,opts.x);
% elseif opts.weighted 
%     [D,Dt] = FD3D_weighted(k,p,q,r,opts.weights); 
else
    [D,Dt] = FD3D_Uniform(k,p,q,r);
end
