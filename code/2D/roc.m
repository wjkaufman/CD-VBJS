%% create ROC curve
% Will Kaufman 2019

clear all;
close all;
addpath('helper_functs');

%% define parameters of run

N = 128;
J = 5;
Jprime = 3;
isChanged = false(1, J); % 1xJ vector that indicates whether measurement
    % is of a "changed" scene
isChanged((Jprime+1):end) = true;

eps = .1;
lam = .25;
order = 2;
willDisp = true;

% ROC curve generation
iter = 2^2; % number of iterations to perform for the ROC curve
numDetect = 2^4; % number of pixels in changed region to sample
numFA = 2^4; % number of pixels in unchanged region to sample
T = 64; % number of threshold values to evaluate at (points along curve)
thresh = reshape(linspace(0, 1, T), [1,1,T]);
thresh = repmat(thresh, N, N, 1);
funct = 'hill';

% spatial oversampling ratio
%(will use os^2 spatial values to inform every frequency value)
os = 2^4;
% level of underdeterminedness (fraction of data kept for analysis)
underdetRatio = 0.95;
std_noise = 1;

%% problem setup

[x,y,f,Y,SNR, f_jump, f_meas, f_VBJS_wl1, changeRegion] = make_data(N, J, ...
    Jprime, funct, order, os, underdetRatio, std_noise, willDisp);


%% GLRT CD

[change, isArrival] = GLRT2D(x, y, isChanged, f_meas, f_VBJS_wl1, 5, willDisp);

%figure; imagesc(change); colorbar;

%%%%%

% assume it's piecewise constant, so TV is sparse
diffMat = -1 * eye(N);
diffMat((N+1):N+1:end) = 1;
diffMat(end,:) = zeros(1,N);
L = diffMat;

% ROC stuff
% want to see how varying threshold T impacts PD and PFA

pd1 = zeros(T, 1);
pd2 = zeros(T, 1);
pfa1 = zeros(T, 1);
pfa2 = zeros(T, 1);

cumChanged = zeros(N, N, T); % (i,j,t) = fraction of time the pixel (i,j)
            % is marked as a change with threshold t

disp(['Signal to noise ratio is ' num2str(SNR)]);

for i = 1:iter
    disp(['iteration ' num2str(i)]);
    
    % randomly select pixels for true change and no change (false alarm)
    ac = randsample(find(changeRegion), numDetect);
    anc = randsample(find(~changeRegion), numFA);
    
    if i == 1
        willDisp = true;
    else
        willDisp = false;
    end
    
    [x,y,f,Y,SNR, f_jump, f_meas, f_VBJS_wl1, changeRegion] = make_data(N, J, ...
    Jprime, funct, order, os, underdetRatio, std_noise, false);
    % below is what I had before for the 1D case
%     [~, ~, ~, Y] = make_data(ref_func, chg_func, N, K, J, Jprime, ...
%         noise, M, prefix, false);
    
    % IDT I need the following
%     [Ghat] = vbjs_reconstruct(N, K, J, Jprime, x, Y, L, prefix, willDisp);

    [change, isArrival] = GLRT2D(x, y, isChanged, f_meas, f_VBJS_wl1, 5, false);
    
    % isObsChanged is NxNxT logical matrix that records if pixel (i,j) is
    % measured as a change with threshold t
    isObsChanged = repmat(change, 1, 1, T);  % NxNxT matrix
    isObsChanged = isObsChanged > thresh;
    cumChanged = cumChanged + isObsChanged;
    % check if estimated changes match with true changes
    % and see how often estimated changes match with no true changes
    
    pd1 = pd1 + reshape(sum(isObsChanged(ac + N^2 * (0:(T-1))), 1) / ...
                           (numel(ac)), T, 1);
    pd2 = pd2 + reshape(sum(~isObsChanged(anc + N^2 * (0:(T-1))), 1) / ...
                            (numel(anc)), T, 1);
    pfa1 = pfa1 + reshape(sum(isObsChanged(anc + N^2 * (0:(T-1))), 1) / ...
                            (numel(anc)), T,1);
    pfa2 = pfa2 + reshape(sum(~isObsChanged(ac + N^2 * (0:(T-1))), 1) / ...
                            (numel(ac)), T,1);
end

cumChanged = cumChanged / iter;
pd1 = pd1 / iter; pd2 = pd2 / iter;
pfa1 = pfa1 / iter; pfa2 = pfa2 / iter;


csvwrite('roc_data.csv', [pfa1, pd1, pfa2, pd2])

% then plot stuff
figure; plot(pfa1, pd1, '-*', [0 1], [0 1], 'k-.');
title('Receiver operator curve (detect = change)');
xlabel('PFA'); ylabel('PD');
set(gcf, 'PaperPosition', [0 0 7 5]);
set(gcf, 'PaperSize', [7 5]);

figure; plot(pfa2, pd2, '-*', [0 1], [0 1], 'k-.');
title('Receiver operator curve (detect = no change)');
xlabel('PFA'); ylabel('PD');
set(gcf, 'PaperPosition', [0 0 7 5]);
set(gcf, 'PaperSize', [7 5]);
