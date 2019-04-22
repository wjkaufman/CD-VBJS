%% create ROC curve
% Will Kaufman 2019

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
disp = false;

% ROC curve generation
iter = 8; % number of iterations to perform for the ROC curve
numDetect = 5; % number of pixels in changed region to sample
numFA = 5; % number of pixels in unchanged region to sample
T = 100; % number of threshold values to evaluate at (points along curve)
thresh = reshape(linspace(0, 1, T), [1,1,T]);
thresh = repmat(thresh, N, N, 1);
funct = 'hill';

os = 2^4; % spatial oversampling ratio 
        %(will use os^2 spatial values to inform every frequency value)
std_noise = 0.5;

%% problem setup

[x,y,f,Y,SNR, f_jump, f_meas, f_VBJS_wl1, changeRegion] = make_data(N, J, ...
    Jprime, funct, order, os, std_noise, disp);


%% GLRT CD

change = GLRT2D(x, y, isChanged, f_meas, f_VBJS_wl1, 5, disp);

%figure; imagesc(change); colorbar;

%%%%%

% assume it's piecewise constant, so TV is sparse
diffMat = -1 * eye(N);
diffMat((N+1):N+1:end) = 1;
diffMat(end,:) = zeros(1,N);
L = diffMat;

% ROC stuff
% want to see how varying threshold T impacts PD and PFA

pd = zeros(T, 1);
pfa = zeros(T, 1);

cumChanged = zeros(N, N, T); % (i,j,t) = fraction of time the pixel (i,j)
            % is marked as a change with threshold t

% OLD CODE (TODO: remove)
% % define where actual changes and actual no changes occured (try multiple)
% ac1 = false(N, T);
% ac1(abs(x) > .2 & abs(x) < .23, :) = true; % change occurs in [.2,.25]
% anc1 = false(N, T);
% anc1(abs(x) > .4, :) = true;
% 
% % define another one, harder
% ac2 = false(N,T);
% ac2(abs(x) > .23 & abs(x) < .25, :) = true;
% anc2 = false(N,T);
% anc2(abs(x) > .25 & abs(x) < .3, :) = true;

disp(['Signal to noise ratio is ' num2str(SNR)]);

for i = 1:iter
    disp(['iteration ' num2str(i)]);
    
    % randomly select pixels for true change and no change (false alarm)
    ac = randsample(find(changeRegion), numDetect);
    anc = randsample(find(~changeRegion), numFA);
    
    if i == 1
        disp = true;
    else
        disp = false;
    end
    
    [x,y,f,Y,SNR, f_jump, f_meas, f_VBJS_wl1, changeRegion] = make_data(N, J, ...
    Jprime, funct, order, os, std_noise, disp);
    % below is what I had before for the 1D case
%     [~, ~, ~, Y] = make_data(ref_func, chg_func, N, K, J, Jprime, ...
%         noise, M, prefix, false);
    
    % IDT I need the following
%     [Ghat] = vbjs_reconstruct(N, K, J, Jprime, x, Y, L, prefix, disp);

    change = GLRT2D(x, y, isChanged, f_meas, f_VBJS_wl1, 5, disp);
    
    isChanged = repmat(change, 1, 1, T);  % NxNxT matrix
    isChanged = isChanged > thresh;
    cumChanged = cumChanged + isChanged;
    % check if estimated changes match with true changes
    % and see how often estimated changes match with no true changes
    pd = pd + reshape(sum(isChanged(mod(ac-1, N)+1, ceil(ac/N), :), ...
                          [1,2]) / numel(ac),...
                      T, 1);
    pfa = pfa + reshape(sum(isChanged(mod(anc-1, N)+1, ceil(anc/N), :), ...
                            [1,2]) / numel(anc), ...
                        T,1);
end

cumChanged = cumChanged / iter;
pd = pd / iter; pfa = pfa / iter;

% then plot stuff
figure; plot(pfa, pd, '-*', [0 1], [0 1], 'k-.');
title('Receiver operator curve');
xlabel('PFA'); ylabel('PD');
set(gcf, 'PaperPosition', [0 0 7 5]);
set(gcf, 'PaperSize', [7 5]);

% figure; imagesc(cumChanged'); title('Not sure what to call this');
% % TODO change the axis numbers to correspond to x/thresh vals
% xlabel('Domain'); ylabel('Threshold values');
% set(gcf, 'PaperPosition', [0 0 7 5]);
% set(gcf, 'PaperSize', [7 5]);
% print([prefix sprintf('idk-N_%d-K_%d-J_%d', N, K, J)], '-dpdf');