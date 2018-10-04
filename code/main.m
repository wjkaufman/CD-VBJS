%% main script for MATH 76 research project
% Will Kaufman 2018
%
%close all; clear all; % apparently MATLAB thinks this is inefficient...

prefix = '../graphics/ex1_';
printGraphs = true;

x0 = 1;
a = 1/2;
ref_func = @(k) toy_func(x0, a, k);
x1 = -.25;
b = 1/4;
chg_func = @(k) toy_func(x1, b, k);

N = 1000;
K = 50;
J = 10;
Jprime = 5;
noise = 1e-2;

[x, SNR, changed, Y] = make_data(ref_func, chg_func, N, K, J, Jprime, ...
    noise, prefix, printGraphs);
xChanged = abs(x) <= b; % logical, is there change? Depends on toy function

% assume it's piecewise constant, so TV is sparse
diffMat = -1 * eye(N);
diffMat((N+1):N+1:end) = 1;
diffMat(end,:) = zeros(1,N);
L = diffMat;

[Ghat] = vbjs_reconstruct(N, K, J, Jprime, x, Y, L, prefix, printGraphs);

[change] = glrt(N, K, J, Jprime, x, changed, Y, Ghat, 5, prefix, printGraphs);