function [change] = GLRT2D(x, y, changed, f_meas, f_VBJS, ...
    nbhdSize)
% do GLRT based on hypothesis testing, assume Gaussian likelihoods
% x and y are coordinates of pixel locations
% f_meas is individually-reconstructed data (both ref and changed)
% f_VBJS is reconstruction of data using VBJS (best estimate of scene)

% get zero mean for GLRT assumption
f_norm = f_meas - f_VBJS;

% if printGraphs
%     figure; plot(x, f_norm, '*'); title('f_norm vs. x');
%     set(gcf, 'PaperPosition', [0 0 7 5]);
%     set(gcf, 'PaperSize', [7 5]);
%     print([prefix sprintf('f_norm-N_%d-K_%d-J_%d', N, K, J)], '-dpdf');
% end

% matrix that records changed regions
change = zeros(length(x), length(y));

J = length(changed); % total number of measurements made
Jp = sum(~changed); % total number of reference measurements

% TODO: implement using blockproc instead?
for i = 1:(length(x)-nbhdSize + 1)
    for j = 1:(length(y)-nbhdSize + 1)
        indX = i:(i+nbhdSize - 1); indY = j:(j+nbhdSize - 1);

        numerator = abs(1/Jp*sum(f_norm(indX, indY,~changed).^2,'all'))^Jp * ...
                    abs(1/(J-Jp)*sum(f_norm(indX, indY, changed).^2,'all'))^(J-Jp);
        denominator = abs(1/J*(sum(f_norm(indX,indY,~changed).^2,'all') + ... % ref images
                               sum(f_norm(indX,indY, changed).^2,'all')))^J; % changed images
        
        % TODO: think about different hypothesis for multiple changed
        % images: see if all changed images are same, or if day1 ? day2
        % etc.
        
        change(i+floor(nbhdSize/2),j+floor(nbhdSize/2)) = ...
            1-(numerator/denominator);
    end
end

% if printGraphs
%     figure; plot(x, f_meas); hold on; title('Y vs. x, with CD');
%     plot(x, change, 'k-', ...
%         'LineWidth', 4);
%     hold off;
%     set(gcf, 'PaperPosition', [0 0 7 5]);
%     set(gcf, 'PaperSize', [7 5]);
%     print([prefix sprintf('Y_CD-N_%d-K_%d-J_%d', N, K, J)], '-dpdf');
% end

end