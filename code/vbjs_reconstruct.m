function [Ghat] = vbjs_reconstruct(N, K, J, Jprime, x, Y, L, ...
    prefix, printGraphs)
% do VBJS reconstruction, stuff and things

yArray1 = repmat(Y, 1, 1, J);
yArray2 = permute(repmat(Y, 1, 1, J), [1 3 2]);
D = sum((yArray1 - yArray2).^2).^(1/2);
D = reshape(D, J, J);

if printGraphs
    figure; imagesc(D); title('Distances between MMVs'); xlabel('Measurement number');
    ylabel('Measurement number'); colorbar;
    set(gcf, 'PaperPosition', [0 0 7 5]);
    set(gcf, 'PaperSize', [7 5]);
    print([prefix sprintf('D-N_%d-K_%d-J_%d', N, K, J)], '-dpdf');
end

% pick best vector for reconstructing reference image

Dref = D(1:Jprime, 1:Jprime);
Dref(1:Jprime+1:end) = Inf;
optimalRefInd = mod(find(Dref == min(min(Dref)), 1, 'first')-1, Jprime)+1;

% VBJS reconstruction

% calculate sparse domain for Y
P = L*Y;
if printGraphs
    figure; plot(x, P); title('P vs. x (sparse domain)');
    set(gcf, 'PaperPosition', [0 0 7 5]);
    set(gcf, 'PaperSize', [7 5]);
    print([prefix sprintf('P-N_%d-K_%d-J_%d', N, K, J)], '-dpdf');
end

% calculate variance across measurements for each spatial location

v = var(P, 0, 2);
if printGraphs
    figure; plot(x, v); title('variance vs. x');
%     set(gcf, 'PaperPosition', [0 0 7 5]);
%     set(gcf, 'PaperSize', [7 5]);
%     print([prefix sprintf('var-N_%d-K_%d-J_%d', N, K, J)], '-dpdf');
end

% calculate weights for weighted l1 norm

PNorm = abs(P)./max(abs(P), [], 1);
C = 1/J * sum(sum(PNorm));

% TODO confirm that this threshold value is good
thresh = 25/N;
indices = 1/J * sum(PNorm, 2) > thresh;
weights = C * (1 - v/max(v));
weights(indices) = 1/C^2 * weights(indices);
if printGraphs
    plot(x, weights); title('weights vs. x');
%     set(gcf, 'PaperPosition', [0 0 7 5]);
%     set(gcf, 'PaperSize', [7 5]);
%     print([prefix sprintf('weights-N_%d-K_%d-J_%d', N, K, J)], '-dpdf');
end

% solve optimization problem for Ghat

Ghat = zeros(N, 1);

l = 1;
cvx_begin
    variable Ghat(N,1)
    minimize(norm(Ghat-Y(:,optimalRefInd), 2) + l*norm(diag(weights) * L * Ghat, 1))
cvx_end

if printGraphs
    figure; plot(x, Ghat); title('Ghat vs. x');
    set(gcf, 'PaperPosition', [0 0 7 5]);
    set(gcf, 'PaperSize', [7 5]);
    print([prefix sprintf('Ghat-N_%d-K_%d-J_%d', N, K, J)], '-dpdf');
end

end