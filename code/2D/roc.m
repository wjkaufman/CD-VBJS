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
iter = 100; % number of iterations to perform for the ROC curve
numDetect = 5; % number of pixels in changed region to sample
numFA = 5; % number of pixels in unchanged region to sample
T = 100; % number of threshold values to evaluate at (points along curve)
thresh = repmat(linspace(0, 1, T), N, 1);

funct = 'hill';

os = 2^4; % spatial oversampling ratio 
        %(will use os^2 spatial values to inform every frequency value)
std_noise = 0.55;

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

cumChanged = zeros(N, T);

% define where actual changes and actual no changes occured (try multiple)
ac1 = false(N, T);
ac1(abs(x) > .2 & abs(x) < .23, :) = true; % change occurs in [.2,.25]
anc1 = false(N, T);
anc1(abs(x) > .4, :) = true;

% define another one, harder
ac2 = false(N,T);
ac2(abs(x) > .23 & abs(x) < .25, :) = true;
anc2 = false(N,T);
anc2(abs(x) > .25 & abs(x) < .3, :) = true;

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
    %
    % TODO as of 2019-04-17 16:26:10
    % Keep on working here, actually sample the change/no change regions
    % and store the relevant values...
    % Then make pretty plots!
    % also try to remember what cumChanged is (I don't remember...)
    isChanged = repmat(change', 1, T);  % N-by-T matrix
    isChanged = isChanged > thresh;
    cumChanged = cumChanged + isChanged;
    % check if estimated changes match with true changes
    % and see how often estimated changes match with no true changes
    pd1 = pd1 + (sum(isChanged .* ac1, 1)/sum(ac1(:,1)))';
    pfa1 = pfa1 + (sum(isChanged .* anc1, 1)/sum(anc1(:,1)))';
    % and do it for the other change/no change regions
    pd2 = pd2 + (sum(isChanged .* ac2, 1)/sum(ac2(:,1)))';
    pfa2 = pfa2 + (sum(isChanged .* anc2, 1)/sum(anc2(:,1)))';
end

cumChanged = cumChanged / iter;
pd1 = pd1 / iter; pfa1 = pfa1 / iter;
pd2 = pd2 / iter; pfa2 = pfa2 / iter;

% then plot stuff
figure; plot(pfa1, pd1, '-*', ...
    pfa2, pd2, '--+', [0 1], [0 1], 'k-.');
title('Receiver operator curve');
xlabel('PFA'); ylabel('PD');
set(gcf, 'PaperPosition', [0 0 7 5]);
set(gcf, 'PaperSize', [7 5]);
print([prefix sprintf('roc-N_%d-K_%d-J_%d', N, K, J)], '-dpdf');

% figure; imagesc(cumChanged'); title('Not sure what to call this');
% % TODO change the axis numbers to correspond to x/thresh vals
% xlabel('Domain'); ylabel('Threshold values');
% set(gcf, 'PaperPosition', [0 0 7 5]);
% set(gcf, 'PaperSize', [7 5]);
% print([prefix sprintf('idk-N_%d-K_%d-J_%d', N, K, J)], '-dpdf');