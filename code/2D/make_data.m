function [x, y, f, Y, SNR, ...
    f_jump, f_meas, f_VBJS_wl1, changeRegion] = make_data(N, J, Jprime, ...
                                    funct, order, os, std_noise, disp)
% returns noisy data with changes, and optionally prints graphs to files
%
%%% inputs %%%
% funct: string that determines function type
% os: oversampling ratio to determine Fourier coefficients
% std_noise: standard deviation of noise in frequency domain
% disp: boolean (to display plots)
%
%%% outputs %%%
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
changeMagnitude = 5;
[X_os,Y_os] = meshgrid(x_os,y_os);
f_chgd = changeMagnitude*(X >= u & X <= (u+du) & ...
                          Y >= v & Y <= (v+dv));
f_chgd_os = changeMagnitude*(X_os >= u & X_os <= (u+du) & ...
                             Y_os >= v & Y_os <= (v+dv));
F_CHGD_os(:,:,(Jprime+1):J) = repmat(f_chgd_os, 1, 1, J-Jprime);
% make change in the first Jprime (so that the "change" is a removal)
% F_CHGD_os(:,:,1:Jprime) = repmat(f_chgd_os, 1, 1, Jprime);

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

% forward operator
% A is Fourier operator, AH is adjoint (= inverse)
A =  @(u) reshape(fft2(reshape(u, N, N)) / sqrt(numel(u)), N^2, 1);
AH = @(u) reshape(ifft2(reshape(u, N, N)) *sqrt(numel(u)), N^2, 1);
% separate operator for oversampling in spatial domain
A_os =  @(u) reshape(fft2(reshape(u, os*N, os*N)) / sqrt(os^2*numel(u)), (os*N)^2, 1);

Y = zeros(N,N,J);
f_meas = zeros(N,N,J);
f_jump = zeros(N,N,J, 2); % approximation to jump function
            % calculated by concentration factors
            % (x,y,measurement, x-dir)
k_vals = [0:N/2 -N/2+1:-1];
l_vals = k_vals'; % might be different if image is not square
k_mat = repmat(k_vals, N, 1);
l_mat = repmat(l_vals, 1, N);
SNR = zeros(1,J);
if disp
    figure; plot(x,f(:,N/2),'k--','linewidth',1.25); hold on;
end

% observe J measurements
for ii = 1:J
    sprintf('on iter %d', ii)
    tmp = f_os+F_CHGD_os(:,:,ii); % add change to data
    Y_os = reshape(A_os(tmp), os*N, os*N) + noise_os(:,:,ii);
    SNR(ii) = snr(Y_os-noise_os(:,:,ii), noise_os(:,:,ii));
    % downsample the Fourier coefficients to only get low-frequency
    % info (what we'd normally get if we measured according to the
    % normal spatial gridpoints
    Y(:,:,ii) = Y_os([1:floor(N/2), ceil(N*(os-.5)+1):os*N], ...
        [1:floor(N/2), ceil(N*(os-.5)+1):os*N]);
    
    f_star = real(reshape(AH(Y(:,:,ii)), N, N)); % do inverse Fourier
    f_meas(:,:,ii) = f_star;
    
    if disp
        plot(x,f_meas(N/2,:,ii),'linewidth',1.25);
    end
    
    % and jump function calculation
    conc_factor_orders = [1, 4, 16];
    jump_x_mat = zeros(N,N, length(conc_factor_orders));
    jump_y_mat = zeros(N,N,length(conc_factor_orders));
    for jj = 1:numel(conc_factor_orders)
        jump_x_mat(:,:,jj) = real(reshape(AH(...
            conc_factor(k_mat, conc_factor_orders(jj)).*Y(:,:,ii)), N, N));
        jump_y_mat(:,:,jj) = real(reshape(AH(...
            conc_factor(l_mat, conc_factor_orders(jj)).*Y(:,:,ii)), N, N));
    end
    f_jump(:,:,ii,1) = minmod(jump_x_mat);
    f_jump(:,:,ii,2) = minmod(jump_y_mat);
end

SNR = mean(SNR, 'all');

if disp
    h = xlabel('$x$');
    xlim([min(x) max(x)]);
    set(h,'interpreter','latex','fontsize',18);
    h = ylabel('$f(x,0)$');
    set(h,'interpreter','latex','fontsize',18);
    set(gca,'fontname','times','fontsize',16);
end

if disp
    % plot jump function in x direction
    figure; colormap gray;
    imagesc(x,y,f_jump(:,:,J,1));
    colorbar; axis xy image;
    title('Jump function approx. in x direction');
    h = xlabel('$x$');
    xlim([min(x) max(x)]);
    set(h,'interpreter','latex','fontsize',18);
    h = ylabel('$y$');
    set(h,'interpreter','latex','fontsize',18);
    set(gca,'fontname','times','fontsize',16);
    
    % jump function in y direction
    figure; colormap gray;
    imagesc(x,y,f_jump(:,:,J,2));
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
% only want to get best _reference_ image, which would be first one
% [~,meas_mat] = get_VWJSdata(reshape(f_meas(:,:,1:(J)),N^2,J));
% [j_star,meas_mat] = get_VWJSdata(...
%     reshape(f_meas(:,:,1:(Jprime)),N^2,Jprime));
% data_js = Y(:,:,j_star);

% manually set optimal data vector for CD
j_star = 1;
data_js = Y(:,:,j_star);

%% variance and weights

% variance
v = var(f_jump,1,3);

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

% take weights to be minimum o W_x, W_y
W = min(W, [], 3);

if disp
    figure; imagesc(x,y,W);
    title('min weight from both');
    axis xy image; colorbar;
    xticks([]);
    yticks([]);
end
%% reconstructions

% VBJS wl1 optimization parameters
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

if disp
    figure;
    colormap gray
    imagesc(x,y,real(f_VBJS_wl1),dyn_range);
    axis xy image; colorbar;
    title('ADMM reconstruction')
    xticks([]);
    yticks([]);
end

end
