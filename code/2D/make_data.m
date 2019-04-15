function [x, SNR, changed, Y] = make_data(ref_func, chg_func, N, J, ...
    Jprime, noise, M, prefix, disp)
% returns noisy data with changes, and optionally prints graphs to files
%
% ref_func: function that returns reference image
% chg_func: function that returns changed image
% N: # samples in spatial domain to reconstruct
% K: # Fourier coefficients gathered
% J: # MMVs gathered
% Jprime: # reference images
% noise: if you don't know what this is, that's too bad
% M: row selector matrix for Fourier coefficients

% select which vectors are "changed" images
changed = false(1,J);
changed(end-Jprime+1:end) = true; % last three are changed images

x = linspace(-1,1,N); % domain of function [-1,1]
fk = repmat(ref_func(0:K)', 1, J); % K-by-J, fourier coefficients
% make change for changed images, simple change
fk(:, changed) = fk(:, changed) + chg_func(0:K)';
% and add noise
signal = mean(fk, 2);
fk = fk + normrnd(0, noise, K+1, J) + 1i * normrnd(0, noise, K+1, J);
SNR = snr(signal, normrnd(0, noise, K+1, 1));
% remove Fourier coefficients according to row selector matrix M
fk = M * fk;

if disp
    figure; plot(real(fk)); title(['Real part of Fourier coefficients' ...
        sprintf(' (SNR=%d)', SNR)]);
    set(gcf, 'PaperPosition', [0 0 7 5]);
    set(gcf, 'PaperSize', [7 5]);
    print([prefix sprintf('fourier_real-N_%d-K_%d-J_%d', N, K, J)], '-dpdf');
end

% individual recovery using partial sum
Fbasis = exp(1i * repmat((0:K)', 1, N) * pi .* x);

% TODO change this assumption, let Y be complex (consider phase)
Y = real(Fbasis' * fk);

if disp
    figure; plot(x, Y); title('Recovered y values');
    set(gcf, 'PaperPosition', [0 0 7 5]);
    set(gcf, 'PaperSize', [7 5]);
    print([prefix sprintf('yhat-N_%d-K_%d-J_%d', N, K, J)], '-dpdf');
end
end