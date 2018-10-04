function [change] = glrt(N, K, J, Jprime, x, changed, Y, Ghat, ...
    nbhdSize, prefix, printGraphs)
% do GLRT based on hypothesis testing, assume Gaussian likelihoods
% Y is individually-reconstructed data (both ref and changed)
% Ghat is reconstruction of data using VBJS (best estimate of scene)

% get zero mean for GLRT assumption
Gnorm = Y - Ghat;

if printGraphs
    figure; plot(x, Gnorm, '*'); title('Gnorm vs. x');
    set(gcf, 'PaperPosition', [0 0 7 5]);
    set(gcf, 'PaperSize', [7 5]);
    print([prefix sprintf('Gnorm-N_%d-K_%d-J_%d', N, K, J)], '-dpdf');
end

change = zeros(size(x));

for i = 1:(length(x)-nbhdSize + 1)
    ind = i:(i+nbhdSize - 1);
    
    numerator = abs(1/Jprime*sum(sum(Gnorm(ind,~changed).^2)))^Jprime * ...
                abs(1/(J-Jprime)*sum(sum(Gnorm(ind, changed).^2)))^(J-Jprime);
    denominator = abs(1/J*(sum(sum(Gnorm(ind,~changed).^2)) + ... % ref images
                           sum(sum(Gnorm(ind,changed).^2))))^J; % changed images
    
    % TODO
    % continue thinking about change statistic, what's appropriate and
    % what's not?
    
    change(i+floor(nbhdSize/2)) = 1- (numerator/denominator);
end

if printGraphs
    figure; plot(x, Y); hold on; title('Y vs. x, with CD');
    plot(x, change, 'k-', ...
        'LineWidth', 4);
    hold off;
    set(gcf, 'PaperPosition', [0 0 7 5]);
    set(gcf, 'PaperSize', [7 5]);
    print([prefix sprintf('Y_CD-N_%d-K_%d-J_%d', N, K, J)], '-dpdf');
end

end