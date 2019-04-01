%% concentration factor based VBJS

clear all
close all

%run('/home/scarnta/Documents/MATLAB/CVX/cvx_startup.m')
addpath(genpath('helper_funct'));

% parameters
N = 128; % resolution 
funct = 'sawtooth_mpi_pi'; % which function to reconstruct
J = 10; % number of measurements
Jprime = 5; % number of reference measurements
dist = 'l2'; % distance function to use 
disp = 1;  % for displaying figures
PA_order = 2; % order of PA transform 
SNR = 20; % signal to noise ratio 

%% Generate data

% Fourier coefficients
k = [0:N/2-1 N/2 -N/2+1:-1]';

% function and grid
[x,f_true] = get_funct(N,funct,k);

% change in domain
f_chg = zeros(size(f_true));
f_chg(floor(N/5):floor(2*N/5)) = mean(abs(f_true));

% forward model
A = dftmtx(N);

% noise
s = mean(abs(A*f_true));
mu = s/SNR;
noise = mu*(randn(N,J) + 1i*randn(N,J));

% initialize
Y = zeros(N,J);
f_tild = zeros(N,J);
sig = zeros(N,J);
order_cf = 2;

if disp
    figure(1); hold on;
    legendtext1 = cell(1,J);
    figure(2); hold on;
    legendtext2 = cell(1,J);
end

for ii = 1:J
    
    if ii <= Jprime
        Y(:,ii) = A*f_true + noise(:,ii);
    else
        Y(:,ii) = A*(f_true + f_chg) + noise(:,ii);
    end
    [f_tild(:,ii),sig(:,ii)] = dual_edge_fhat_multiorder_old(order_cf(end),N, Y(:,ii));
    
    if disp
        figure(1);
        h1 = plot(x,f_tild(:,ii));
        legendtext1{ii} = ['Meas = ',num2str(ii)];
        
        figure(2);
        h2 = plot(x,sig(:,ii));
        legendtext2{ii} = ['Order = ',num2str(order_cf(end))];
    end
    order_cf = [order_cf; order_cf + 2];
    
end
if disp
    figure(1); legend(legendtext1)
    figure(2); legend(legendtext2)
end
%% VBJS

f_jump_ref = f_tild(:,1:Jprime);

[w,v] = get_VBJSweights_CF(f_jump_ref);

if strcmp(dist,'l2')
    [j_star,meas_mat] = get_VBJSdata_l2(f_jump_ref);
elseif strcmp(dist,'emd')
    [j_star,meas_mat] = get_VBJSdata_emd(f_jump_ref);
elseif strcmp(dist,'mhd')
    [j_star,meas_mat] = get_VBJSdata_mh(f_jump_ref);
end

data_js = Y(:,j_star);
if disp
    figure; imagesc(meas_mat)
    figure; plot(x,w,'-x');
end
%% final reconstructions

% PA operators
PA = PA_Operator_1D( N,PA_order );
W = diag(w);

% l1 reg
cvx_begin quiet
clear cvx
variable f_star_l1(N,1)
minimize( norm(W*PA*f_star_l1,1) + norm(A*f_star_l1 - data_js,2));
cvx_end

% l2 reg
f_star_l2 = (A'*A+PA'*(W'*W)*PA)\(A'*data_js);

l2_error= norm(f_star_l2-f_true,2)/norm(f_true,2);
l1_error = norm(f_star_l1-f_true,2)/norm(f_true,2);

if disp
    leg1 = sprintf('CF VBJS $\\ell_1$: Error = %7.4f\n',l1_error);
    leg2 = sprintf('CF VBJS $\\ell_2$: Error = %7.4f\n',l2_error);
    
    figure; hold on;
    plot(x,f_true,'--','linewidth',1.5)
    plot(x,f_star_l1,':','linewidth',1.5);
    plot(x,f_star_l2,'-.','linewidth',1.5);
    L = legend('True',leg1,leg2);
    set(L,'interpreter','latex','fontsize',14);
    h = xlabel('$x$');
    set(h,'interpreter','latex','fontsize',18);
    h = ylabel('$f(x)$');
    set(h,'interpreter','latex','fontsize',18);
    set(gca,'fontname','times','fontsize',14);
    %         ylim([-1.1,1.75])
    xlim([min(x),max(x)])
end

%% CD

% individually reconstruct data using inverse Fourier transform
f_meas = zeros(N, J);
for ii = 1:J
    f_meas(:,ii) = real(ifft(Y(:,ii))); % I think... TODO check
end
if disp
    figure; plot(f_meas(:,[1, 6]));
    title('individual reconstructions (measurement 1 and 6)');
end

[gamma] = glrt_1d(J, Jprime, x, f_meas, f_star_l1, 5);

if disp
    figure; plot(gamma); title('change statistic');
end