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
J = 5;
Jprime = 3;
eps = .1; 
lam = .25; 
order = 2;

funct = 'hill';

os = 2^4;
std_noise = 0.55;

%% problem setup 

[x,y,f,Y,SNR, f_jump, f_meas, f_VBJS_wl1, changeRegion] = make_data(N, J, ...
    Jprime, funct, order, os, std_noise, false);


%% GLRT CD

% make changed vector that records which measurements are "changed"
changed = false(1, J);
changed((Jprime+1):end) = true;

change = GLRT2D(x, y, changed, f_meas, f_VBJS_wl1, 5, 1);

figure; imagesc(change); colorbar;
