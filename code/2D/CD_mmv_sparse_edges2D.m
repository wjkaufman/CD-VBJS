%%%
%
% Modifications (adding CD algorithm at the end) by Will Kaufman (2018-2019)
%
% code and helper functions from Theresa Scarnati (2018)
%
%%%

clear all;
close all;
addpath('helper_functs');

%% parameters

N = 128;
num_meas = 5;
num_chgd = 2;
eps = .1; 
lam = .25; 
order = 2;

funct = 'hill';

%% problem setup 

% true image
f = get_img(funct,N);
dyn_range = [min(min(f)),max(max(f))];
std_noise = .55;%std(f(:))*1e-1;

x = linspace(-1,1,N);
y = linspace(-1,1,N);

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

% noise
noise = std_noise*randn(N,N,num_meas);% + 1i*std_noise*randn(N,N,M);

% changed data, goes at the front of the list of measurements
F_CHGD = zeros(N,N,num_meas);
% make consistent change for last num_chgd measurements
u = -.75; du = .1;
v = -.1; dv = .1;
[X,Y] = meshgrid(x,y);
f_chgd = 5*(X >= u & X <= (u+du) & Y >= v & Y <= (v+dv)); % I think fine
F_CHGD(:,:,(num_meas-num_chgd+1):num_meas) = repmat(f_chgd, 1, 1, num_chgd);
figure(20); colormap gray;
imagesc(x,y,f_chgd+f,dyn_range);
colorbar; axis xy image;
h = xlabel('$x$');
xlim([min(x) max(x)]);
set(h,'interpreter','latex','fontsize',18);
h = ylabel('$y$');
set(h,'interpreter','latex','fontsize',18);
set(gca,'fontname','times','fontsize',16);

% add false information too (one step at a time though)
% TODO _eventually_ add false information

% forward operator
% A is Fourier operator, AH is adjoint (= inverse)
A =  @(u) reshape(fft2(reshape(u, N, N)) / sqrt(numel(u)), N^2, 1); 
AH = @(u) reshape(ifft2(reshape(u, N, N)) *sqrt(numel(u)), N^2, 1);

% PA operator 
PA = PA_Operator_1D(N,order);

% noisy data and measurements 
% opts.outer_iter = 50; 
% opts.inner_iter = 20; 
% opts.scale_b = false; 
% opts.scale_A = false; 
% opts.weighted = false; 
% opts.data_mlp = true; 

Y = zeros(N,N,num_meas);
f_meas = zeros(N,N,num_meas);
f_jump = zeros(N,N,num_meas, 2); % approximation to jump function
            % calculated by concentration factors
            % (x,y,measurement, x-dir)
k_vals = [0:N/2 -N/2+1:-1];
l_vals = k_vals'; % might be different if image is not square
k_mat = repmat(k_vals, N, 1);
l_mat = repmat(l_vals, 1, N);
PAf_meas = zeros(N,N,num_meas);
PAf_meas_vec = zeros(N^2,num_meas);
figure; plot(x,f(:,N/2),'k--','linewidth',1.25); hold on;
for ii = 1:num_meas
    sprintf('on iter %d', ii)
    tmp = f+F_CHGD(:,:,ii); % add change to data
    Y(:,:,ii) = reshape(A(tmp), N, N) + noise(:,:,ii);
    
    f_star = real(reshape(AH(Y(:,:,ii)), N, N)); % do inverse Fourier sum
    f_meas(:,:,ii) = f_star;
    
    % and sparse domain calculation
    jump_x = real(reshape(AH(...
        conc_factor(k_mat).*Y(:,:,ii)), N, N));
    jump_y = real(reshape(AH(...
        conc_factor(l_mat).*Y(:,:,ii)), N, N));
    f_jump(:,:,ii,1) = jump_x;
    f_jump(:,:,ii,2) = jump_y;
    
    % old code (TODO still need?)
    PAf_meas(:,:,ii) = real(PA*f_star + f_star*PA);
    PAf_meas_vec(:,ii) = col(PAf_meas(:,:,ii));
    plot(x,f_meas(N/2,:,ii),'linewidth',1.25);
end
h = xlabel('$x$');
xlim([min(x) max(x)]);
set(h,'interpreter','latex','fontsize',18);
h = ylabel('$f(x,0)$');
set(h,'interpreter','latex','fontsize',18);
set(gca,'fontname','times','fontsize',16);

% TODO 
% filter jump functions by comparing sign across all
% measurement vectors: if same sign -> keep value
% if different sign -> set jump to 0

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

%% optimal data vector
% only want to get best _reference_ image, so just look through that
% TODO: or do I just want to manually set j_star = 1? A priori reason
% that this should be our reference (we're calling this t=0)...
% I'm manually setting j_star, so replaced it in line below with `~`
% [~,meas_mat] = get_VWJSdata(reshape(f_meas(:,:,1:(num_meas)),N^2,num_meas));
% [j_star,meas_mat] = get_VWJSdata(...
%     reshape(f_meas(:,:,1:(num_meas - num_chgd)),N^2,num_meas-num_chgd));
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

% plot variance in x direction
figure; imagesc(x,y,v(:,:,1)); title('variance in x direction');
colorbar; axis xy image;
h = xlabel('$x$');
xlim([min(x) max(x)]);
set(h,'interpreter','latex','fontsize',18);
h = ylabel('$y$');
set(h,'interpreter','latex','fontsize',18);
set(gca,'fontname','times','fontsize',16);

% weights
wx = get_VWJSweights(reshape(f_jump(:,:,:,1), N^2, num_meas),.15);
wy = get_VWJSweights(reshape(f_jump(:,:,:,2), N^2, num_meas),.15);
W = zeros(N,N, 2);
W(:,:,1) = reshape(wx, N, N);
W(:,:,2) = reshape(wy, N, N);

figure; imagesc(x,y,W(:,:,1));
axis xy image; colorbar;
xticks([]); 
yticks([]);

figure; imagesc(x,y,W(:,:,2));
axis xy image; colorbar;
xticks([]); 
yticks([]);

% for now, just take the minimum of the weights (not sure
%   how to do ADMM with two terms to minimize)
% TODO figure this out!

W = min(W, [], 3);
figure; imagesc(x,y,W);
title('min weight from both');
axis xy image; colorbar;
xticks([]); 
yticks([]);

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

% TODO this reconstruction is not correct
[f_VBJS_wl1,out_wl1] = ADMM2(A,AH,data_js,[N,N],opts_wl1);

figure;
colormap gray
imagesc(x,y,real(f_VBJS_wl1),dyn_range);
axis xy image; colorbar;
title('ADMM reconstruction')
xticks([]); 
yticks([]);
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
%% GLRT CD

% make changed vector that records which measurements are "changed"
changed = false(1, num_meas);
changed((num_meas - num_chgd+1):end) = true;

change = GLRT2D(x, y, changed, f_meas, f_VBJS_wl1, 5, 1);

figure; imagesc(change); colorbar;
