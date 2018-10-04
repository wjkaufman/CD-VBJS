%% VBJS reconstruction

prefix = '../graphics/vbjs01-';

N = 101; % # samples in spatial domain to reconstruct
K = 50; % # Fourier coefficients gathered
J = 5; % # MMVs gathered
Jprime = 5;
noise = 1e-2;

% make row selector matrix
M = eye(K+1);

% make sparsifying transform
diffMat = -1 * eye(N);
diffMat((N+1):N+1:end) = 1;
diffMat(end,:) = zeros(1,N);
L = diffMat;

x0 = 1;
a = 1/2;
ref_func = @(k) toy_func(x0, a, k);

[x, SNR, changed, Y] = make_data(ref_func, ref_func, N, K, J, Jprime, ...
        noise, M, prefix, true);
[Ghat] = vbjs_reconstruct(N, K, J, Jprime, x, Y, L, prefix, true);