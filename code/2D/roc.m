%% create ROC curve
% Will Kaufman 2019

% define parameters of run

prefix = '../../graphics/undet_01_1-';

N = 128; % number of spatial points to reconstruct
K = N; % number of Fourier coefficients observed (normally K = N)
J = 5; % total number of measurements made on scene
Jprime = 3; % total number of rereference measurements (no change)

noise = 1e-3; % TODO change to SNR

% TODO modify the existing code below to work for 2D
% I'll need to refactor the code to easily generate/oversample spatial
%   domain, then find Fourier coefficients, etc. etc.

x0 = 1;
a = 1/2;
ref_func = @(k) toy_func(x0, a, k);
x1 = -.25;
b = 1/4;
chg_func = @(k) toy_func(x1, b, k);

[x, SNR, changed, ~] = make_data(ref_func, chg_func, N, K, J, Jprime, ...
        noise, M, prefix, true);
xChanged = abs(x) <= b; % logical, is there change? Depends on toy function

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

% assume it's piecewise constant, so TV is sparse
diffMat = -1 * eye(N);
diffMat((N+1):N+1:end) = 1;
diffMat(end,:) = zeros(1,N);
L = diffMat;

% ROC stuff
% want to see how varying threshold T impacts PD and PFA
iter = 100;
T = 200; % number of threshold values to evaluate at
thresh = repmat(linspace(0, 1, T), N, 1);

pd1 = zeros(T, 1);
pfa1 = zeros(T, 1);
pd2 = zeros(T, 1);
pfa2 = zeros(T, 1);
cumChanged = zeros(N, T);

disp(['Signal to noise ratio is ' num2str(SNR)]);

for i = 1:iter
    disp(['iteration ' num2str(i)]);
    if i == 1
        printGraphs = true;
    else
        printGraphs = false;
    end
%     M = diag(rand(1,K+1) < .8); % randomly do it every time
    [~, ~, ~, Y] = make_data(ref_func, chg_func, N, K, J, Jprime, ...
        noise, M, prefix, false);
    
    [Ghat] = vbjs_reconstruct(N, K, J, Jprime, x, Y, L, prefix, printGraphs);

    [change] = glrt(N, K, J, Jprime, x, changed, Y, Ghat, 3, prefix, printGraphs);

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