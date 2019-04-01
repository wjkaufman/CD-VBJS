function [ f,out ] = grad_descent_mmv(A,AH,y,W,N,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs:
%   A    matrix or function handle representing forward operator (p x q)
%   AH   transpose of forward operator (q x p)
%   y    data matrix (p x r)
%   W    weighting matrix (q x r)
%   N    contains size information of solution [q,r]
%   opts options
% outputs:
%   f    final, optimized solution (q x r)
%   out  output values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initializations

global P p q r


if ~isa(A,'function_handle')
    A_mat = A;
    [A,AH] = get_handles(A,N);
end

[p,r] = size(y);
q = N(1);

if isfield(opts,'f_init')
    f = opts.f_init;
else
    f = AH(y);
end

f = f(:);
y = y(:);

tol = opts.tol; %1e-6;
max_it = opts.max_it; %100;
max_bt = opts.max_bt; %25;
rel_chg = 1e3;
cnt = 1;
delta = opts.delta; %1e-4;
rho = opts.rho; %.4;
tau = opts.tau; % 1
lam = opts.lam;

P = PA_Operator_1D(q,opts.order);

if ~isfield(opts,'isreal')
    opts.isreal = 1;
    theta = ones(q,q);
elseif opts.isreal == 0 
    theta = exp(-1i*angle(reshape(f,q,q)));
elseif opts.isreal == 1
    theta = ones(q,q); 
end

if ~isfield(opts,'dyn_range')
    dyn_range = [min(min(f)),max(max(f))];
else
    dyn_range = opts.dyn_range;
end

% get obj function
J = get_J(A,f,y,W,lam,theta);

% initial step length
out.f = f(:);
out.rel_chg = [];
out.alpha = [];
out.obj = J;

if r == 1
    f_goal = (lam*P'*W'*W*P+A_mat'*A_mat)\(A_mat'*y);
end


%% gradient descent
flg = 1;
while and((rel_chg > tol), (cnt < max_it))
    
    dJ = get_dJ(A,AH,f,y,W,lam,theta);
    
    % set step length
    if flg
        % first iteration
        alpha = tau;
        flg = 0;
    else
        s = f - fp;
        u = dJ - dJp;
        
        % bb step length
        alpha = abs((s'*s)/(s'*u));
    end
    
    f_tmp = f-alpha*dJ;
    J = get_J(A,f,y,W,lam,theta);
    J_l = get_J(A,f_tmp,y,W,lam,theta);
    J_r = J - delta*alpha*(dJ(:)')*dJ(:);
    
    iter_bt = 1;
    % backtrack while Armijo cond not satisfied
    while and((J_l > J_r), (iter_bt < max_bt))
        
        alpha = alpha*rho;
        
        f_tmp = f - alpha*dJ;
        if ~opts.isreal
            theta_tmp = exp(-1i*angle(reshape(f_tmp,q,q)));
        else 
            theta_tmp = theta; 
        end
        J_l = get_J(A,f_tmp,y,W,lam,theta_tmp);
        J_r = J - delta*alpha*(dJ(:)')*dJ(:);
        
        iter_bt = iter_bt + 1;
        
    end
    
    fp = f;
    dJp = dJ;
    
    % update solution and gradient vectors
    f = fp - alpha*dJp;
    
    if ~opts.isreal
        theta = exp(-1i*angle(reshape(f,q,q)));
    end
    
    J = get_J(A,f,y,W,lam,theta);
    

    
    rel_chg = norm(f(:)-fp(:))/norm(fp(:));
    
    out.f = [out.f,f(:)];
    out.alpha = [out.alpha;alpha];
    out.rel_chg = [out.rel_chg;rel_chg];
    out.obj = [out.obj;J];
    
    cnt = cnt + 1;
    
    % plot results
    if opts.disp
        if r == 1
            figure(11);
            if isfield(opts,'f_true')
                plot(opts.f_true,'k'); hold on;
            end
            plot(f_goal,'bo'); hold on;
            plot(f,'rx'); hold off;
            title(strcat('iteration = ',num2str(cnt)));
            if isfield(opts,'f_true')
                legend('true function','true l2 min','alg');
            else
                legend('true l2 min','alg');
            end
            pause(0.1);
        else
            if opts.isreal
                figure(11); imagesc(abs(reshape(f,q,r)),dyn_range);
            else
                figure(11); imagesc(img_dB(reshape(f,q,q)),dyn_range);
            end
            colorbar
            title(strcat('iteration = ',num2str(cnt)));
            pause(0.1)
        end
    end
end

if cnt == max_it
    out.inner_msg = 'max iterations reached';
else
    out.inner_msg = 'relative change below threshold';
end

f = reshape(f,q,r);

end

function J = get_J(A,x,y,W,lam,theta)
global P p q r
% J(x) = 1/2 ||Ax - y||_2^2 + lam/2 ||Px||_2,w^2
% A - p x q
% x - q x r
% y - p x r
% W - q x r

x = double(reshape(x,q,r));
y = double(reshape(y,p,r));

if r == 1
    j1 = .5*norm(W*P*(theta.*x),2)^2;
    j2 = .5*norm(A(x)-y(:),2)^2;
else
    Ax = reshape(A(x),p,r);
    j1 = .5*norm(W.*(P*real(theta.*x)),'fro')^2 +...
        .5*norm(W.*(real(theta.*x)*P),'fro')^2;
    j2 = .5*norm(Ax-y,'fro')^2;
end

J = lam*j1+j2;

end

function G = get_dJ(A,AH,x,y,W,lam,theta)
global P p q r
% d_x J(x) = lam(P'(W.*Px) + (W.*xP)P') + A'(Ax-y)
% A - p x q
% x - q x r
% y - p x r
% W - q x r

x = double(reshape(x,q,r));
y = double(reshape(y,p,r));

if r == 1
%     g1 = lam*P'*W'*W*P*(theta.*x);
    g1 = lam*P'*W'*W*P*(x);
else
    Px = P*real(theta.*x);
    xP = real(theta.*x)*P;
    g1 = lam*(P'*(W.*Px)+(W.*xP)*P');
end


Axy = reshape(A(x),p,r) - y;
g2 = AH(Axy);

G = g1(:) + g2(:);

end