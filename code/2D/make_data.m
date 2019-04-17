function [x, y, f, Y, SNR, ...
    f_jump, f_meas, f_VBJS_wl1, changeRegion] = make_data(N, J, Jprime, funct, order, os, std_noise, disp)
% returns noisy data with changes, and optionally prints graphs to files
%
% funct: string that determines function type
% os: oversampling ratio to determine Fourier coefficients
% std_noise: standard deviation of noise in frequency domain
% disp: boolean (to display plots)
%
% outputs
% x/y: spatial gridpoints on which measurements are made
% f: NxN true image of underlying scene
% Y: NxNxJ Fourier coefficients
% SNR: signal to noise ratio (in db)
% f_jump: NxNxJx2 approximation to jump function, both in x- and y-dir
% f_meas: individual reconstructions using inverse Fourier transform
% f_VBJS_wl1: VBJS reconstruction using l1 regularization
% changeRegion: NxN logical for where the change actually was

% true image
f = get_img(funct,N);
f_os = get_img(funct,os*N);
dyn_range = [min(min(f)),max(max(f))];

x = linspace(-1,1,N);
y = linspace(-1,1,N);
x_os = linspace(-1,1,os*N);
y_os = linspace(-1,1,os*N);

if disp
    figure;
    colormap gray;
    imagesc(x,y,f,dyn_range);
    colorbar; axis xy image;
    h = xlabel('$x$');
    xlim([min(x) max(x)]);
    set(h,'interpreter','latex','fontsize',18);
    h = ylabel('$y$');
    set(h,'interpreter','latex','fontsize',18);
    set(gca,'fontname','times','fontsize',16);
end

% noise
%noise = std_noise*randn(N,N,J);% + 1i*std_noise*randn(N,N,M);
noise_os = std_noise*randn(os*N,os*N,J);

% changed data, goes at the front of the list of measurements
%F_CHGD = zeros(N,N,J);
F_CHGD_os = zeros(os*N,os*N,J);
% make consistent change for last J-Jprime measurements
u = -.75; du = .1;
v = -.1; dv = .1;
[X,Y] = meshgrid(x,y);
% define where the change region is
changeRegion = (X >= u & X <= (u+du) & Y >= v & Y <= (v+dv));
[X_os,Y_os] = meshgrid(x_os,y_os);
f_chgd = 5*(X >= u & X <= (u+du) & Y >= v & Y <= (v+dv)); % I think fine
f_chgd_os = 5*(X_os >= u & X_os <= (u+du) & Y_os >= v & Y_os <= (v+dv)); % I think fine
%F_CHGD(:,:,(Jprime+1):J) = repmat(f_chgd, 1, 1, J-Jprime);
F_CHGD_os(:,:,(Jprime+1):J) = repmat(f_chgd_os, 1, 1, J-Jprime);

if disp
    figure(20); colormap gray;
    imagesc(x,y,f_chgd+f,dyn_range);
    colorbar; axis xy image;
    h = xlabel('$x$');
    xlim([min(x) max(x)]);
    set(h,'interpreter','latex','fontsize',18);
    h = ylabel('$y$');
    set(h,'interpreter','latex','fontsize',18);
    set(gca,'fontname','times','fontsize',16);
end

% add false information too (one step at a time though)
% TODO _eventually_ add false information

% forward operator
% A is Fourier operator, AH is adjoint (= inverse)
A =  @(u) reshape(fft2(reshape(u, N, N)) / sqrt(numel(u)), N^2, 1); 
AH = @(u) reshape(ifft2(reshape(u, N, N)) *sqrt(numel(u)), N^2, 1);
% separate operator for oversampling in spatial domain
A_os =  @(u) reshape(fft2(reshape(u, os*N, os*N)) / sqrt(os^2*numel(u)), (os*N)^2, 1); 

% PA operator 
PA = PA_Operator_1D(N,order);

% noisy data and measurements 
% opts.outer_iter = 50; 
% opts.inner_iter = 20; 
% opts.scale_b = false; 
% opts.scale_A = false; 
% opts.weighted = false; 
% opts.data_mlp = true; 

Y = zeros(N,N,J);
f_meas = zeros(N,N,J);
f_jump = zeros(N,N,J, 2); % approximation to jump function
            % calculated by concentration factors
            % (x,y,measurement, x-dir)
k_vals = [0:N/2 -N/2+1:-1];
l_vals = k_vals'; % might be different if image is not square
k_mat = repmat(k_vals, N, 1);
l_mat = repmat(l_vals, 1, N);
PAf_meas = zeros(N,N,J);
PAf_meas_vec = zeros(N^2,J);
SNR = zeros(1,J);
if disp
    figure; plot(x,f(:,N/2),'k--','linewidth',1.25); hold on;
end
for ii = 1:J
    sprintf('on iter %d', ii)
    tmp = f_os+F_CHGD_os(:,:,ii); % add change to data
    Y_os = reshape(A_os(tmp), os*N, os*N) + noise_os(:,:,ii);
    % TODO how to compute SNR??
    SNR(ii) = 0;% snr(norm(Y_os-noise_os(:,:,ii)), noise_os(:,:,ii));
    % downsample the Fourier coefficients to only get low-frequency
    % info (what we'd normally get if we measured according to the 
    % normal spatial gridpoints
    Y(:,:,ii) = Y_os([1:floor(N/2), ceil(N*(os-.5)+1):os*N], ...
        [1:floor(N/2), ceil(N*(os-.5)+1):os*N]);
    
    f_star = real(reshape(AH(Y(:,:,ii)), N, N)); % do inverse Fourier
    f_meas(:,:,ii) = f_star;
    
    % old code (TODO still need?)
    PAf_meas(:,:,ii) = real(PA*f_star + f_star*PA);
    PAf_meas_vec(:,ii) = col(PAf_meas(:,:,ii));
    if disp
        plot(x,f_meas(N/2,:,ii),'linewidth',1.25);
    end
    
    % and sparse domain calculation
    jump_x = real(reshape(AH(...
        conc_factor(k_mat).*Y(:,:,ii)), N, N));
    jump_y = real(reshape(AH(...
        conc_factor(l_mat).*Y(:,:,ii)), N, N));
    f_jump(:,:,ii,1) = jump_x;
    f_jump(:,:,ii,2) = jump_y;
end
if disp
    h = xlabel('$x$');
    xlim([min(x) max(x)]);
    set(h,'interpreter','latex','fontsize',18);
    h = ylabel('$f(x,0)$');
    set(h,'interpreter','latex','fontsize',18);
    set(gca,'fontname','times','fontsize',16);
end

% TODO
% filter jump functions by comparing sign across all
% measurement vectors: if same sign -> keep value
% if different sign -> set jump to 0

if disp
    % plot jump function
    figure; colormap gray;
    imagesc(x,y,f_jump(:,:,5,1));
    colorbar; axis xy image;
    title('Jump function approx. in x direction');
    h = xlabel('$x$');
    xlim([min(x) max(x)]);
    set(h,'interpreter','latex','fontsize',18);
    h = ylabel('$y$');
    set(h,'interpreter','latex','fontsize',18);
    set(gca,'fontname','times','fontsize',16);

    figure; colormap gray;
    imagesc(x,y,f_jump(:,:,5,2));
    colorbar; axis xy image;
    title('Jump function approx. in y direction');
    h = xlabel('$x$');
    xlim([min(x) max(x)]);
    set(h,'interpreter','latex','fontsize',18);
    h = ylabel('$y$');
    set(h,'interpreter','latex','fontsize',18);
    set(gca,'fontname','times','fontsize',16);
end

%% optimal data vector
% only want to get best _reference_ image, so just look through that
% TODO: or do I just want to manually set j_star = 1? A priori reason
% that this should be our reference (we're calling this t=0)...
% I'm manually setting j_star, so replaced it in line below with `~`
% [~,meas_mat] = get_VWJSdata(reshape(f_meas(:,:,1:(J)),N^2,J));
% [j_star,meas_mat] = get_VWJSdata(...
%     reshape(f_meas(:,:,1:(Jprime)),N^2,Jprime));
% data_js = Y(:,:,j_star);

% manually set optimal data vector for CD
% TODO this should be changed, yes?
j_star = 1;
data_js = Y(:,:,j_star);

% figure; imagesc(meas_mat);
% colorbar;
% h = xlabel('measurement number');
% set(h,'interpreter','latex','fontsize',18);
% h = ylabel('measurement number');
% set(h,'interpreter','latex','fontsize',18);
% set(gca,'fontname','times','fontsize',16);

% figure;
% colormap gray
% imagesc(x,y,real(f_meas(:,:,j_star)),dyn_range);
% axis xy image
% xticks([]); 
% yticks([]);

%% variance and weights

% variance
v = var(f_jump,1,3);
%v = reshape(v,N,N); 

if disp
    % plot variance in x direction
    figure; imagesc(x,y,v(:,:,1)); title('variance in x direction');
    colorbar; axis xy image;
    h = xlabel('$x$');
    xlim([min(x) max(x)]);
    set(h,'interpreter','latex','fontsize',18);
    h = ylabel('$y$');
    set(h,'interpreter','latex','fontsize',18);
    set(gca,'fontname','times','fontsize',16);
end

% weights
wx = get_VWJSweights(reshape(f_jump(:,:,:,1), N^2, J),.15);
wy = get_VWJSweights(reshape(f_jump(:,:,:,2), N^2, J),.15);
W = zeros(N,N, 2);
W(:,:,1) = reshape(wx, N, N);
W(:,:,2) = reshape(wy, N, N);

if disp
    figure; imagesc(x,y,W(:,:,1));
    axis xy image; colorbar;
    xticks([]); 
    yticks([]);

    figure; imagesc(x,y,W(:,:,2));
    axis xy image; colorbar;
    xticks([]); 
    yticks([]);
end

% for now, just take the minimum of the weights (not sure
%   how to do ADMM with two terms to minimize)
% TODO figure this out!

W = min(W, [], 3);
if disp
    figure; imagesc(x,y,W);
    title('min weight from both');
    axis xy image; colorbar;
    xticks([]); 
    yticks([]);
end
%% reconstructions

%W = ones(N,N);

% VBJS wl1
opts_wl1.mu =1 ;
opts_wl1.beta = 1; 
opts_wl1.outer_iter = 50;
opts_wl1.inner_iter = 20;
opts_wl1.scale_b = true;
opts_wl1.scale_A = true;
opts_wl1.weighted = true;
opts_wl1.weights = W;
opts_wl1.data_mlp = true;
opts_wl1.disp = 0;
opts_wl1.order = order;

[f_VBJS_wl1,~] = ADMM2(A,AH,data_js,[N,N],opts_wl1);
% not storing out_wl1 (second output)

if disp
    figure;
    colormap gray
    imagesc(x,y,real(f_VBJS_wl1),dyn_range);
    axis xy image; colorbar;
    title('ADMM reconstruction')
    xticks([]); 
    yticks([]);
end
% %%
% % VBJS wl2
% opts_wl2.max_it = 100; 
% opts_wl2.tol = 1e-10; 
% opts_wl2.max_bt = 25; 
% opts_wl2.delta = 1e-5; 
% opts_wl2.rho = .4; 
% opts_wl2.order = order; 
% opts_wl2.tau = 1; 
% opts_wl2.lam = 1; 
% opts_wl2.disp = 0; 
% opts_wl2.dyn_range = dyn_range; 
% opts_wl2.f_init = A_mat\data_js; 
% 
% [f_VBJS_wl2,out_wl2] = grad_descent_mmv(A,AH,data_js,10*W/max(max(W)),[N,N],opts_wl2); 
% 
% figure; colormap gray
% imagesc(x,y,real(f_VBJS_wl2),dyn_range);
% axis xy image
% xticks([]); 
% yticks([]);

end