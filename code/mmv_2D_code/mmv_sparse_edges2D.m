clear all;
close all;
addpath('helper_functs');

%% parameters

N = 128; 
num_meas = 5; 
num_false = 0; 
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
colormap gray
imagesc(x,y,f,dyn_range);
colorbar; axis xy image
h = xlabel('$x$');
xlim([min(x) max(x)]);
set(h,'interpreter','latex','fontsize',18);
h = ylabel('$y$');
set(h,'interpreter','latex','fontsize',18);
set(gca,'fontname','times','fontsize',16);

% noise
noise = std_noise*randn(N,N,num_meas);% + 1i*std_noise*randn(N,N,M);

% false data
F_FALSE = zeros(N,N,num_meas);
for ii = 1:num_false
    u = -1 + 2*rand(1); % random number between -1 and 1
    [X,Y] = meshgrid(x,y);
    f_false = randi([round(dyn_range(1)), round(dyn_range(2))],1).*(X<=u) +...
        randi([round(dyn_range(1)), round(dyn_range(2))],1).*(Y>u);
    F_FALSE(:,:,ii) = f_false;
    figure(20); imagesc(f_false+f,dyn_range); pause(0.4);
end

% forward operator
A  = randn(N,N);
d = 1./sqrt(sum(A.^2,2));
A = A.*d; A_mat = A; 
[A,AH] = get_handles(A,N);

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
PAf_meas = zeros(N,N,num_meas);
PAf_meas_vec = zeros(N^2,num_meas);
figure; plot(x,f(:,N/2),'k--','linewidth',1.25); hold on;
for ii = 1:num_meas
    tmp = f+F_FALSE(:,:,ii);
    Y(:,:,ii) = reshape(A(tmp),N,N) + noise(:,:,ii);
    
    opts.mu =  1;%randi([1,10],1); % data
    opts.beta = 1;%randi([1,10],1); % l1
    
    [f_star,out] = ADMM2(A,AH,Y(:,:,ii),[N,N],opts); 
    f_meas(:,:,ii) = f_star;
    PAf_meas(:,:,ii) = real(PA*f_star + f_star*PA);
    PAf_meas_vec(:,ii) = col(PAf_meas(:,:,ii));
    plot(x,f_meas(:,N/2,ii),'linewidth',1.25); 
end
h = xlabel('$x$');
xlim([min(x) max(x)]);
set(h,'interpreter','latex','fontsize',18);
h = ylabel('$f(x,0)$');
set(h,'interpreter','latex','fontsize',18);
set(gca,'fontname','times','fontsize',16);

%% optimal data vector
[j_star,meas_mat] = get_VWJSdata(reshape(f_meas,N^2,num_meas));
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
colorbar; axis xy image
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
%%
% VBJS wl2
opts_wl2.max_it = 100; 
opts_wl2.tol = 1e-10; 
opts_wl2.max_bt = 25; 
opts_wl2.delta = 1e-5; 
opts_wl2.rho = .4; 
opts_wl2.order = order; 
opts_wl2.tau = 1; 
opts_wl2.lam = 1; 
opts_wl2.disp = 0; 
opts_wl2.dyn_range = dyn_range; 
opts_wl2.f_init = A_mat\data_js; 

[f_VBJS_wl2,out_wl2] = grad_descent_mmv(A,AH,data_js,10*W/max(max(W)),[N,N],opts_wl2); 

figure; colormap gray
imagesc(x,y,real(f_VBJS_wl2),dyn_range);
axis xy image
xticks([]); 
yticks([]);
%% compare results 

marker = {'k--','-.','-',':'};

figure; plot(x,f(:,64),marker{1},'linewidth',1.5)
hold on; plot(x,real(f_meas(:,64,j_star)),marker{4},'linewidth',1.5)
hold on; plot(x,real(f_VBJS_wl1(64,:)),marker{2},'linewidth',1.5)
hold on; plot(x,real(f_VBJS_wl2(:,64)),marker{3},'linewidth',1.5); 
h = xlabel('$x$');
xlim([min(x) max(x)]);
set(h,'interpreter','latex','fontsize',18);
h = ylabel('$f(x,0)$');
set(h,'interpreter','latex','fontsize',18);
set(gca,'fontname','times','fontsize',16);

L = legend('True Function','SMV','VBJS $\ell_1$','VBJS $\ell_2$');
set(L,'interpreter','latex','fontsize',14,'location','south');
