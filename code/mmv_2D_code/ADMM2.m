function [U,out] = ADMM2(A,AH,b,n,opts)

% Modifications by Theresa Scarnati @AFRL
% 05/24/2018


% Modifications by Toby Sanders @ASU
% School of Math & Stat Sciences
% 08/24/2016


% This code has been modified to solve l1 penalty problems with
% higher order TV operators.  Several small bugs and notation
% changes have been made as well.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Problem Description       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [U, out] = HOTV3D(A,b,n,opts)

% Motivation is to find:

%               min_f { mu/2*||Af - b||_2^2 + ||D^k f||_1 }

% where D^k is kth order finite difference.
% Multiscale finite differences D^k can also be used.
% To see how to modify these settings read the file "check_HOTV_opts.m"

% The problem is modified using variable splitting
% and this algorithm solves:

%      min_{f,w} {mu/2 ||Af - b||_2^2 + beta/2 ||D^k f - w ||_2^2
%               + ||w||_1 - (delta , Af - b ) - (sigma , D^k f - w) }

% delta and sigma are Lagrange multipliers
% Algorithm uses alternating direction minimization over f and w.


% This algorithm was originally authored by Chengbo Li at Rice University
% as a TV solver called TVAL3.
% original code and description can be found here:
% http://www.caam.rice.edu/~optimization/L1/TVAL3/

% Inputs:
%   A: matrix operator as either a matrix or function handle
%   b: data values in vector form
%   n: image/ signal dimensions in vector format
%   opts: structure containing input parameters,
%       see function check_HOTV_opts.m for these


% Outputs:
%   U: reconstructed signal
%   out: output numerics

%%

if numel(n)<3
    n(end+1:3) = 1;
elseif numel(n)>3
    error('n can have at most 3 dimensions');
end
p = n(1); q = n(2); r = n(3);

% get and check opts
opts = check_ADMM_opts(opts);

% mark important constants
tol_inn = opts.tol_inn;
tol_out = opts.tol_out;
k = opts.order;
n = p*q*r;
wrap_shrink = opts.wrap_shrink;
if round(k)~=k
    wrap_shrink = true;
end


% unify implementation of A
if ~isa(A,'function_handle')
    [A,AH] = get_handles(A,p);
end

%check that A* is true adjoint of A
% [flg,~,~] = check_D_Dt(@(u)A(u),@(u)AH(u),[n,1]);
% if ~flg
%     warning('A and A* do not appear consistent');
% end
% clear flg;

% initialize
b = b(:);
b_orig = b; 
% out.data = [b];

% check scaling A
if opts.scale_A
    [A,AH,b] = ScaleA(n,A,AH,b);
end

% check scaling b
scl = 1;
if opts.scale_b
    [b,scl] = Scaleb(b);
end

% check for maximum constraint value
if opts.max_c
    max_v = opts.max_v*scl;
end

% calculate A'*b
Atb = AH(b);

% initialize everything else
global D Dt
[U,mu,beta,muf,betaf,muDbeta,sigma,delta,gL,ind,out] ...
    = get_ADMM(p,q,r,Atb,scl,opts,k,b,wrap_shrink);    % U: p*q

nrmb = norm(b);
Upout = U;

if sum(sum(sum(opts.phase_angles)))==0
    theta = 1;
else
    theta = exp(-1i*opts.phase_angles);
end

Uc = D(U,theta);

%%
% first shrinkage step
if opts.weighted
    W = max(abs(Uc) -  opts.weights(:), 0).*sign(Uc);
else
    W = max(abs(Uc) - 1/beta, 0).*sign(Uc);
end
% reset edge values if not using periodic regularization
if ~wrap_shrink, W(ind)=Uc(ind); end

lam1 = sum(col(abs(W)));

Au = A(U(:));

%%

% gA and gD are the gradients of ||Au-b||^2 and ||Du-w||^2, respectively
% i.e. g = A'(Au-b), gD = D'(Du-w)
[lam2,lam3,lam4,lam5,f,gD,Au,gA] = get_grad(U,Uc,W,...
    lam1,beta,mu,A,AH,b,Atb,sigma,delta,theta,opts);


% compute gradient
g = gD + muDbeta*gA - gL;

% update output vectors
out.f = [out.f; f];
out.lam1 = [out.lam1; lam1]; out.lam2 = [out.lam2; lam2];
out.lam3 = [out.lam3; lam3];out.lam4 = [out.lam4; lam4];
out.lam5 = [out.lam5; lam5];out.mu = [out.mu; mu];
out.sparsity = [];
out.DU = [out.DU;norm(Uc(:),1)];
out.data = [b_orig, b(:)];

rel_chg_out = 0;
for ii = 1:opts.outer_iter
    if opts.disp
        fprintf('    Beginning outer iteration #%d\n',ii);
        fprintf('    mu = %d , beta = %d , order = %g, rel chg =%g\n',mu,beta,k,rel_chg_out);
        fprintf('iter    ||w||_1    ||Du - w||^2  ||Au - b||^2   rel chg\n');
    end
    
    %initialize the constants
    gam = opts.gam; Q = 1; fp = f;
    
    
    for jj = 1:opts.inner_iter
        % compute step length, tau
        if jj~=1
            % BB-like step length
            dgA = gA - gAp;
            dgD = gD - gDp;
            ss = uup'*uup;
            sy = uup'*(dgD + muDbeta*dgA);
            tau = abs(ss/max(sy,eps));
        else
            % do Steepest Descent at the 1st ieration
            gc = D(reshape(g,p,q,r),theta);
            dDd = sum(col(gc.*conj(gc)));
            Ag = A(g);
            tau = abs((g'*g)/(dDd + muDbeta*(Ag')*Ag));
        end
        
        % keep previous values for backtracking & computing next tau
        Up = U; gAp = gA; gDp = gD; Aup = Au;
        Ucp = Uc; %DtsAtdp =  DtsAtd;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ONE-STEP GRADIENT DESCENT %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        U = U(:) - tau*g;
        % projected gradient method for inequality constraints
        if opts.nonneg
            U = max(real(U),0);
        elseif opts.isreal
            U = real(U);
        end
        if opts.max_c
            U = min(U,max_v);
        end
        U = reshape(U,p,q,r);
        Uc = D(U,theta);
        
        [lam2,lam3,lam4,lam5,f,gD,Au,gA] = get_grad(U,Uc,W,...
            lam1,beta,mu,A,AH,b,Atb,sigma,delta,theta,opts);
        
        % Nonmonotone Line Search Back tracking
        % Unew = Up + alpha*(U - Up)
        % f should be decreasing, if not, then the algorithm moves U
        % back in the direction of the previous solution
        alpha = 1;
        du = U - Up;
        if opts.weighted
            const = 1e-5*(g'*g*tau);
        else
            const = 1e-5*beta*(g'*g*tau);
        end
        cnt = 0; flg = true;
        
        while f > fp - alpha*const
            if cnt <5
                if flg
                    dgA = gA - gAp;
                    dgD = gD - gDp;
                    dAu = Au - Aup;
                    dUc = Uc - Ucp;
                    flg = false;
                end
                % shrink alpha
                alpha = alpha*opts.gamma;
                % U is moved back toward Up, in particular:
                % U = alpha*U +(1-alpha)Up;
                % all other values are updated accordingly
                [U,lam2,lam3,lam4,lam5,f,Uc,Au,gA,gD] = back_up(p,q,r,...
                    lam1,alpha,beta,mu,Up,du,gAp,dgA,gDp,dgD,Aup,dAu,W,...
                    Ucp,dUc,b,sigma,delta,opts);
                cnt = cnt + 1;
            else
                
                % shrink gam
                gam = opts.rate_gam*gam;
                
                % give up and take Steepest Descent step
                %                 if (opts.disp > 0) && (mod(jj,opts.disp) == 0)
                %                     disp('    count of back tracking attains 5 ');
                %                 end
                
                % compute step length, tau
                gc = D(reshape(g,p,q,r),theta);
                dDd = sum(col(gc.*conj(gc)));
                Ag = A(g);
                tau = abs((g'*g)/(dDd + muDbeta*(Ag')*Ag));
                %update
                U = Up(:) - tau*g;
                % projected gradient method for inequality constraints
                if opts.nonneg
                    U = max(real(U),0);
                elseif opts.isreal
                    U = real(U);
                end
                
                U = reshape(U,p,q,r);
                Uc = D(U,theta);
                % shrinkage
                if opts.weighted
                    Ucbar = Uc - sigma;
                    W = max(abs(Ucbar) -  opts.weights(:), 0).*sign(Ucbar);
                else
                    Ucbar = Uc - sigma/beta;
                    W = max(abs(Ucbar) - 1/beta, 0).*sign(Ucbar);
                end
                
                % reset edge values if not using periodic regularization
                if ~wrap_shrink, W(ind)=Uc(ind); end
                
                lam1 = sum(col(abs(W)));
                [lam2,lam3,lam4,lam5,f,gD,Au,gA] = get_grad(U,Uc,W,...
                    lam1,beta,mu,A,AH,b,Atb,sigma,delta,theta,opts);
                alpha = 0; % remark the failure of back tracking
                break;
            end
            
            if sum(sum(sum(opts.phase_angles)))~=0
                theta = exp(-1i*angle(reshape(U,p,q,r)));
            end
            
        end
        
        
        
        % if back tracking is successful, then recompute
        if alpha ~= 0
            if opts.weighted
                Ucbar = Uc - sigma;
                W = max(abs(Ucbar) -  opts.weights(:), 0).*sign(Ucbar);
            else
                Ucbar = Uc - sigma/beta;
                W = max(abs(Ucbar) - 1/beta, 0).*sign(Ucbar);
            end
            % reset edge values if not using periodic regularization
            if ~wrap_shrink, W(ind)=Uc(ind); end
            % update parameters related to Wx, Wy
            [lam1,lam2,lam4,f,gD] = update_W(beta,...
                W,Uc,sigma,lam1,lam2,lam4,f,theta,opts);
        end
        
        % update reference value
        Qp = Q; Q = gam*Qp + 1; fp = (gam*Qp*fp + f)/Q;
        uup = U - Up; uup = uup(:);           % uup: pqr
        rel_chg_inn = norm(uup)/norm(Up(:));
        
        
        
        out.f = [out.f; f]; out.C = [out.C; fp]; out.cnt = [out.cnt;cnt];
        out.lam1 = [out.lam1; lam1]; out.lam2 = [out.lam2; lam2]; out.lam3 = [out.lam3; lam3];
        out.lam4 = [out.lam4; lam4]; out.lam5 = [out.lam5; lam5];
        out.tau = [out.tau; tau]; out.alpha = [out.alpha; alpha];out.mu = [out.mu; mu];
        out.rel_chg_inn = [out.rel_chg_inn;rel_chg_inn];
        out.rel_lam2 = [out.rel_lam2;sqrt(lam2)/norm(W(:))];
        out.DU = [out.DU; norm(Uc(:),1)];
        if opts.store_soln
            out.Uall(:,:,jj+(ii-1)*opts.inner_iter) = U;
        end        
        
        if (opts.disp > 0) && (mod(ii,opts.disp) == 0)
            prnt_format = '%3.0f %10.5g %12.5g %13.5g %10.5f\n';
            fprintf(prnt_format, jj,lam1,lam2,lam3,rel_chg_inn);%,out.DU(end));
        end
        
        
        % add new gradient terms
        g = gD + muDbeta*gA - gL;
        
        % move to next outer iteration and update multipliers if relative
        % change is less than tolerance
        if (rel_chg_inn < tol_inn), break; end;
        
        sparsity = norm(Uc -W,1); 
        out.sparsity = [out.sparsity; sparsity];
    end
    % end of inner loop
    
    Au = A(U(:));
    
    % update Atb
    Atb = AH(b);
    out.data = [out.data, b(:)];
    
    rel_chg_out = norm(U(:)-Upout(:))/norm(Upout(:));
    out.rel_chg_out = [out.rel_chg_out; rel_chg_out];
    Upout = U;
    
    % stop if already reached optimal solution
    if rel_chg_out < tol_out || sqrt(lam3(end))/nrmb<opts.min_l2_error
        break;
    end
    
    % update multipliers
    deltap = delta;
    lam5p = lam5;
    [sigma,delta,lam4,lam5] = update_mlp(beta,mu, ...
        W,Uc,Au,b,sigma,delta,opts);
    if ~opts.data_mlp
        delta = deltap; lam5 = lam5p;
    end
    
    
    % update penality parameters for continuation scheme
    %beta0 = beta;

    if opts.weighted
        beta = 1;
        mu = 1;
        muDbeta = 1;
        % update function value, gradient, and relavent constant
        f = lam1 + .5*lam2 + .5*lam3 - lam4 - lam5;
        %gL = -(beta0/beta)*g;     % DtsAtd should be divided by new beta
        gL = (Dt(sigma,theta) + AH(delta));
    else
        beta = min(betaf, beta*opts.rate_ctn);
        mu = min(muf, mu*opts.rate_ctn);
        muDbeta = mu/beta;
        % update function value, gradient, and relavent constant
        f = lam1 + beta/2*lam2 + mu/2*lam3 - lam4 - lam5;
        %gL = -(beta0/beta)*g;     % DtsAtd should be divided by new beta
        gL = 1/beta*(Dt(sigma,theta) + AH(delta));
    end
    % gradient, divided by beta
    g = gD + muDbeta*gA - gL;
    
%     figure(11); imagesc(img_dB(reshape(U,p,q,r)),[-64,0]); pause(0.5);
    
end

out.total_iter = numel(out.f)-1;
out.final_error = norm(A(U(:))-b)/nrmb;
out.final_wl1 = lam1(end);
out.final_Du_w = lam2(end);
out.rel_error = sqrt(out.lam3)/nrmb;
out. mu = mu; out.beta = beta;

if opts.disp
    if out.rel_error(end) < opts.min_l2_error
        fprintf('\nREACHED OPTIMAL L2 ERROR!!!\n\n');
    end
    
    if opts.disp_conv
        final_disp_new(out,opts);
    end
    
end

% rescale U
U = U/scl;

U = reshape(U,p,q,r);


function [lam2,lam3,lam4,lam5,f,gD,Au,gA] = get_grad(U,Uc,W,...
    lam1,beta,mu,A,AH,b,Atb,sigma,delta,theta,opts)
global Dt

Au = A(U(:));

% gA = A'(Au-b)
gA = AH(Au) - Atb;

% lam2, ||Du - w||^2
V = Uc - W;
lam2 = sum(col(V.*conj(V)));

% gD = D'(Du-w)
gD = Dt(V,theta);

% lam3, ||Au - b||^2
Aub = Au-b;
lam3 = Aub'*Aub;%norm(Aub)^2;

%lam4
lam4 = sum(col(sigma.*V));

%lam5
lam5 = delta'*Aub;

if opts.weighted
    % f
    f = lam1 + .5*lam2 + .5*lam3 - lam4 - lam5;
else
    % f
    f = lam1 + beta/2*lam2 + mu/2*lam3 - lam4 - lam5;
end


function [U,lam2,lam3,lam4,lam5,f,Uc,Au,gA,gD] = back_up(p,q,r,lam1,...
    alpha,beta,mu,Up,du,gAp,dgA,gDp,dgD,Aup,dAu,W,Ucp,dUc,...
    b,sigma,delta,opts)

gA = gAp + alpha*dgA;
gD = gDp + alpha*dgD;
U = Up + alpha*reshape(du,p,q,r);
Au = Aup + alpha*dAu;
Uc = Ucp + alpha*dUc;

V = Uc - W;


lam2 = sum(col(V.*conj(V)));
Aub  = Au-b;
lam3 = norm(Aub)^2;
lam4 = sum(col(sigma.*V));
lam5 = delta'*Aub;

if opts.weighted 
    f = lam1 + .5*lam2 + .5*lam3 - lam4 - lam5;
else
    f = lam1 + beta/2*lam2 + mu/2*lam3 - lam4 - lam5;
end

function [lam1,lam2,lam4,f,gD] = update_W(beta,...
    W,Uc,sigma,lam1,lam2,lam4,f,theta,opts)
global Dt

% update parameters because W was updated
if opts.weighted 
    tmpf = f -lam1 - .5*lam2 + lam4;
else
    tmpf = f -lam1 - beta/2*lam2 + lam4;
end
lam1 = sum(col(abs(W)));
V = Uc - W;

gD = Dt(V,theta);
lam2 = sum(col(V.*conj(V)));
lam4 = sum(col(sigma.*V));
if opts.weighted
    f = tmpf +lam1 + .5*lam2 - lam4;
else
    f = tmpf +lam1 + beta/2*lam2 - lam4;
end

function [sigma,delta,lam4,lam5] = update_mlp(beta,mu, ...
    W,Uc,Au,b,sigma,delta,opts)


V = Uc - W;
Aub = Au-b;
if opts.weighted 
    delta = delta - Aub;
    sigma = sigma - V;
else
    sigma = sigma - beta*V;
    delta = delta - mu*Aub;
end
%tmpf = f + lam4 + lam5;
lam4 = sum(col(sigma.*V));
lam5 = delta'*Aub;
%f = tmpf - lam4 - lam5;




