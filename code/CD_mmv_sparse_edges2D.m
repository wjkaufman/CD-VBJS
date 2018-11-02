%%%
%
% Modifications (adding CD algorithm at the end) by Will Kaufman (2018)
%
% code and helper functions from Theresa Scarnati (2018)
%
%%%

%clear all;
%close all;
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

% forward operator
A = @(u) fft2(u) / sqrt(numel(u)); % A is Fourier operator
AH = @(u) sqrt(numel(u)) * ifft2(u); % because apparently adjoint of Fourier operator
                % is inverse

% PA operator 
PA = PA_Operator_1D(N,order);

% noisy data and measurements 
opts.outer_iter = 50; 
opts.inner_iter = 20; 
opts.scale_b = false; 
opts.scale_A = false; 
opts.weighted = false; 
opts.data_mlp = true; 

Y = zeros(N,N,num_meas);
f_meas = zeros(N,N,num_meas);
f_jump = zeros(N,N,num_meas); % approximation to jump function
            % calculated by concentration factors
PAf_meas = zeros(N,N,num_meas);
PAf_meas_vec = zeros(N^2,num_meas);
figure; plot(x,f(:,N/2),'k--','linewidth',1.25); hold on;
for ii = 1:num_meas
    sprintf('on iter %d', ii)
    tmp = f+F_CHGD(:,:,ii); % add change to data
    Y(:,:,ii) = A(tmp) + noise(:,:,ii);
    
    % old ADMM reconstructions, not doing that bc I want
    % a linear transformation (to preserve noise)
    %opts.mu =  1;%randi([1,10],1); % data
    %opts.beta = 1;%randi([1,10],1); % l1
    %[f_star,out] = ADMM2(A,AH,Y(:,:,ii),[N,N],opts);
    
    f_star = real(AH(Y(:,:,ii))); % do simple inverse Fourier sum
    f_meas(:,:,ii) = f_star;
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

%% optimal data vector
% only want to get best _reference_ image, so just look through that
% TODO: or do I just want to manually set j_star = 1? A priori reason
% that this should be our reference (we're calling this t=0)...
% I'm manually setting j_star, so replaced it in line below with `~`
[~,meas_mat] = get_VWJSdata(reshape(f_meas(:,:,1:(num_meas)),N^2,num_meas));
% [j_star,meas_mat] = get_VWJSdata(...
%     reshape(f_meas(:,:,1:(num_meas - num_chgd)),N^2,num_meas-num_chgd));
% data_js = Y(:,:,j_star);

% manually set optimal data vector for CD
j_star = 1;
data_js = Y(:,:,j_star);

figure; imagesc(meas_mat);
colorbar;
h = xlabel('measurement number');
set(h,'interpreter','latex','fontsize',18);
h = ylabel('measurement number');
set(h,'interpreter','latex','fontsize',18);
set(gca,'fontname','times','fontsize',16);

figure;
colormap gray
imagesc(x,y,real(f_meas(:,:,j_star)),dyn_range);
axis xy image
xticks([]); 
yticks([]);

%% variance and weights

% variance
v = var(PAf_meas_vec,1,2);
v = reshape(v,N,N); 

figure; imagesc(x,y,v);
colorbar; axis xy image;
h = xlabel('$x$');
xlim([min(x) max(x)]);
set(h,'interpreter','latex','fontsize',18);
h = ylabel('$y$');
set(h,'interpreter','latex','fontsize',18);
set(gca,'fontname','times','fontsize',16);

% weights
w = get_VWJSweights(PAf_meas_vec,.15);
W = reshape(w,N,N);

figure; imagesc(x,y,W);
axis xy image
xticks([]); 
yticks([]);

%% reconstructions

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

[f_VBJS_wl1,out_wl1] = ADMM2(A,AH,data_js,[N,N],opts_wl1);

figure;
colormap gray
imagesc(x,y,real(f_VBJS_wl1),dyn_range);
axis xy image
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

change = GLRT2D(x, y, changed, f_meas, f_VBJS_wl1, 5);

figure; imagesc(change); colorbar;
