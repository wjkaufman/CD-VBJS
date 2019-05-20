function [x, y, f, Y, SNR, ...
    f_jump, f_meas, f_VBJS_wl1, changeRegion] = make_data(N, J, Jprime, ...
                                    funct, order, os, underdetRatio, std_noise, willDisp)
% returns noisy data with changes, and optionally prints graphs to files
%
%%% inputs %%%
% N: number of spatial values per dimension (so 2D image is N^2 pixels)
% J: total number of measurements
% J': total number of reference measurements (no change)
% funct: string that determines function type
% os: oversampling ratio to determine Fourier coefficients
% undetRatio: level of underdeterminedness (0.8 = 80% of Fourier
%   coefficients are kept)
% std_noise: standard deviation of noise in frequency domain
% disp: boolean to display plots
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
numTargets = 64;
f = get_img_pttarget(N, numTargets, 0);
f_os = get_img_pttarget(os*N, numTargets, 0);
% TODO fix this, this will make two different changed images...
f_chg_os = get_img_pttarget(os*N, numTargets, .25);
f_chg = get_img_pttarget(N, numTargets, .25); % change 25% of pt targets
changeRegion = f_chg ~= 0; % changed region, simple for point target code
dyn_range = [min(min(f)),max(max(f))];

x = linspace(-1,1,N);
y = linspace(-1,1,N);
x_os = linspace(-1,1,os*N);
y_os = linspace(-1,1,os*N);

if willDisp
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

if willDisp
    figure(20); colormap gray;
    imagesc(x,y,f_chg,dyn_range);
    colorbar; axis xy image;
    h = xlabel('$x$');
    xlim([min(x) max(x)]);
    set(h,'interpreter','latex','fontsize',18);
    h = ylabel('$y$');
    set(h,'interpreter','latex','fontsize',18);
    set(gca,'fontname','times','fontsize',16);
end

% forward, inverse operators
% Fourier is Fourier operator, invFourier is adjoint (= inverse)
Fourier =  @(u) reshape(fft2(reshape(u, N, N)) / sqrt(numel(u)), N^2, 1);
invFourier = @(u) reshape(ifft2(reshape(u, N, N)) *sqrt(numel(u)), N^2, 1);
% separate operator for oversampling in spatial domain
Fourier_os =  @(u) reshape(fft2(reshape(u, os*N, os*N)) / sqrt(os^2*numel(u)), (os*N)^2, 1);

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
if willDisp
    figure; plot(x,f(:,N/2),'k--','linewidth',1.25); hold on;
end

% observe J measurements
for ii = 1:J
    disp(['on measurement ', num2str(ii)]);
    if ii <= Jprime % reference image
        tmp = f_os;
    else % changed image
        tmp = f_chg_os;
    end
    Y_os = reshape(Fourier_os(tmp), os*N, os*N) + noise_os(:,:,ii);
    SNR(ii) = snr(Y_os-noise_os(:,:,ii), noise_os(:,:,ii));
    % downsample the Fourier coefficients to only get low-frequency
    % info (what we'd normally get if we measured according to the
    % normal spatial gridpoints
    Y(:,:,ii) = Y_os([1:floor(N/2), ceil(N*(os-.5)+1):os*N], ...
        [1:floor(N/2), ceil(N*(os-.5)+1):os*N]);
    % remove Fourier coefficients randomly to create underdetermined system
    underdetMat = zeros(N, N);
    underdetMat(randsample(N^2, round(N^2 * underdetRatio))) = 1;
    Y(:,:,ii) = underdetMat .* Y(:,:,ii);
    
    f_star = real(reshape(invFourier(Y(:,:,ii)), N, N)); % do inverse Fourier
    f_meas(:,:,ii) = f_star;
    
    if willDisp
        plot(x,f_meas(N/2,:,ii),'linewidth',1.25);
    end
    
    % and jump function calculation
    conc_factor_orders = [1, 4, 16];
    jump_x_mat = zeros(N,N, length(conc_factor_orders));
    jump_y_mat = zeros(N,N,length(conc_factor_orders));
    for jj = 1:numel(conc_factor_orders)
        jump_x_mat(:,:,jj) = real(reshape(invFourier(...
            conc_factor(k_mat, conc_factor_orders(jj)).*Y(:,:,ii)), N, N));
        jump_y_mat(:,:,jj) = real(reshape(invFourier(...
            conc_factor(l_mat, conc_factor_orders(jj)).*Y(:,:,ii)), N, N));
    end
    f_jump(:,:,ii,1) = minmod(jump_x_mat);
    f_jump(:,:,ii,2) = minmod(jump_y_mat);
end

SNR = mean(SNR, 'all');

if willDisp
    h = xlabel('$x$');
    xlim([min(x) max(x)]);
    set(h,'interpreter','latex','fontsize',18);
    h = ylabel('$f(x,0)$');
    set(h,'interpreter','latex','fontsize',18);
    set(gca,'fontname','times','fontsize',16);
end

if willDisp
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

if willDisp
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

if willDisp
    figure; imagesc(x,y,W(:,:,1));
    axis xy image; colorbar; title('weights for x direction');
    xticks([]);
    yticks([]);

    figure; imagesc(x,y,W(:,:,2)); 
    axis xy image; colorbar; title('weights for x direction');
    xticks([]);
    yticks([]);
end

% take weights to be minimum o W_x, W_y
W = min(W, [], 3);

if willDisp
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

[f_VBJS_wl1,~] = ADMM2(Fourier,invFourier,data_js,[N,N],opts_wl1);

if willDisp
    figure;
    colormap gray
    imagesc(x,y,real(f_VBJS_wl1),dyn_range);
    axis xy image; colorbar;
    title('ADMM reconstruction')
    xticks([]);
    yticks([]);
end

end
