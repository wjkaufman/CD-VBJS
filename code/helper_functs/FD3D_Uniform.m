function [D,Dt] = FD3D_Uniform(k,p,q,r)
% this one should be used to form the PA operator on uniform points
% when U is not periodic.

% Written by Toby Sanders @ASU
% School of Math & Stat Sciences
% 08/24/2016

% finite difference operators for polynomial annihilation
% k is the order of the PA transform
D = @(U,theta)D_Forward(U,theta,k,p,q,r);
Dt = @(Uc,theta)D_Adjoint(Uc,theta,k,p,q,r);

% high order finite differences
function [dU] = D_Forward(U,theta,k,p,q,r)
U = reshape(U,p,q,r);
U = U.*theta;
if k~=0
    dU = zeros(p,q,r,3);
    if k<=q
        dU(:,:,:,1) = diff([U,U(:,1:k,:)],k,2);
    end
    if k<=p
        dU(:,:,:,2) = diff([U;U(1:k,:,:)],k,1);
    end
    if k<=r
        dU(:,:,:,3) = diff(cat(3,U,U(:,:,1:k)),k,3);
    end
    dU = reshape(dU,p*q*r,3);
else
    %standard l1 minimization of order 0
    dU = U(:);
end
% 
% dU = [zeros(1,3); dU(2:end-1,:); zeros(1,3)];
dU = dU*2^(1-k);  % normalization

%transpose FD
function U = D_Adjoint(dU,theta,k,p,q,r)
if k~=0
    U = zeros(p,q,r);
    dU = reshape(dU,p,q,r,3);
    if k<=q
        U = U + (-1)^k*diff([dU(:,end-k+1:end,:,1),dU(:,:,:,1)],k,2);
    end    
    if k<=p
        U = U + (-1)^k*diff([dU(end-k+1:end,:,:,2);dU(:,:,:,2)],k,1);
    end
    if k<=r
        U = U + (-1)^k*diff(cat(3,dU(:,:,end-k+1:end,3),dU(:,:,:,3)),k,3);
    end
    U = U(:);
else
    %standard l1 minimization
    U = dU(:);
end
U = U.*conj(theta(:));
U = U*2^(1-k);  % normalization
    
    
    
    
