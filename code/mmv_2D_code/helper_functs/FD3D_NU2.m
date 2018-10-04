function [D,Dt] = FD3D_NU2(k,p,q,r,x)
% thing you are trying to reconstruct is size p x q x r.
% k is the order of HOTV.
% U is the think you want to reconstruct
% Uc is the conjugate of U

% finite difference operators for polynomial annihilation
% k is the order of the PA transform
D = @(U,theta)D_Forward(U,theta,k,p,q,r,x);
Dt = @(Uc,theta)D_Adjoint(Uc,theta,k,p,q,r,x);


% high order finite differences
function [dU] = D_Forward(U,theta,k,p,q,r,xw)

U = reshape(U,p,q);
U = U.*theta;
dU = zeros(p,q,r,3);

x = col(xw);
N = length(x);

if mod(k,2) == 0
    flg = 1;
    m2 = k/2;
else
    flg = 0;
    m12 = (k+1)/2;
end


P = zeros(N,N);


for j = 1:N
    
    if flg
        ind = j-m2:j+m2;
    else
        ind = j-m12:j+m12-1;
    end
    
    if max(ind) > N
        
        I = find(ind<=N);
        ind = ind(I);
        m = length(ind);
        ind_mid = ceil(m/2);
        x_tmp = x(ind);
        wi = zeros(1,m);
        
        for i = 1:m
            poly = x_tmp(i)-x_tmp;
            poly(i) = [];
            wi(i) = prod(poly);
        end
        
        C_spa = factorial(m-1)./wi ;
        C_spa = C_spa/sum(C_spa(1:ind_mid));
        P(j,ind) = C_spa;
        
    elseif min(ind) <= 0
        
        I = find(ind>0);
        ind = ind(I);
        m = length(ind);
        ind_mid = ceil(m/2);
        x_tmp = x(ind);
        wi = zeros(1,m);
        
        for i = 1:m
            poly = x_tmp(i)-x_tmp;
            poly(i) = [];
            wi(i) = prod(poly);
        end
        
        C_spa = factorial(m-1)./wi ;
        C_spa = C_spa/sum(C_spa(1:ind_mid));
        P(j,ind) = C_spa;
        
    else
        
        m = length(ind);
        ind_mid = ceil(m/2);
        x_tmp = x(ind);
        wi = zeros(1,m);
        
        for i = 1:m
            poly = x_tmp(i)-x_tmp;
            poly(i) = [];
            wi(i) = prod(poly);
        end
        
        C_spa = factorial(m)./wi ;
        C_spa = C_spa/sum(C_spa(1:ind_mid));
        P(j,ind) = C_spa;
        
    end
    
    if length(ind) == 1
        P(j,ind) = 0;
    end
    
end

P = sparse(P);%*2^(1-m);

if q ==1
    dU(:,:,:,2) = P*U;
else
    
    dU(:,:,:,2) = P*U+U*P';
end

dU = reshape(dU,p*q*r,3);
dU = dU*2^(1-k);  % normalization

%transpose FD
function U = D_Adjoint(dU,theta,k,p,q,r,xw)

U = zeros(p,q,r);
dU = reshape(dU,p,q,r,3);

x = col(xw);
N = length(x);

if mod(k,2) == 0
    flg = 1;
    m2 = k/2;
else
    flg = 0;
    m12 = (k+1)/2;
end

P = zeros(N,N);

for j = 1:N
    
    if flg
        ind = j-m2:j+m2;
    else
        ind = j-m12:j+m12-1;
    end
    
    if max(ind) > N
        
        I = find(ind<=N);
        ind = ind(I);
        m = length(ind);
        ind_mid = ceil(m/2);
        x_tmp = x(ind);
        wi = zeros(1,m);
        
        for i = 1:m
            poly = x_tmp(i)-x_tmp;
            poly(i) = [];
            wi(i) = prod(poly);
        end
        
        C_spa = factorial(m-1)./wi ;
        C_spa = C_spa/sum(C_spa(1:ind_mid));
        P(j,ind) = C_spa;
        
    elseif min(ind) <= 0
        
        I = find(ind>0);
        ind = ind(I);
        m = length(ind);
        ind_mid = ceil(m/2);
        x_tmp = x(ind);
        wi = zeros(1,m);
        
        for i = 1:m
            poly = x_tmp(i)-x_tmp;
            poly(i) = [];
            wi(i) = prod(poly);
        end
        
        C_spa = factorial(m-1)./wi ;
        C_spa = C_spa/sum(C_spa(1:ind_mid));
        P(j,ind) = C_spa;
        
    else
        
        m = length(ind);
        ind_mid = ceil(m/2);
        x_tmp = x(ind);
        wi = zeros(1,m);
        
        for i = 1:m
            poly = x_tmp(i)-x_tmp;
            poly(i) = [];
            wi(i) = prod(poly);
        end
        
        C_spa = factorial(m)./wi ;
        C_spa = C_spa/sum(C_spa(1:ind_mid));
        P(j,ind) = C_spa;
        
    end
    
    if length(ind) == 1
        P(j,ind) = 0;
    end
    
end

P = sparse(P);%*2^(1-m);

U = P'*dU(:,2);

U = U(:).*conj(theta(:));
U = U*2^(1-k);  % normalization

