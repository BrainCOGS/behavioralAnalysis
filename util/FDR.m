function [isSig, alpha_correct] = FDR(pvals, alpha)

% [isSig, alpha_correct] = FDR(pvals, alpha)
%
% performs a false discovery rate correction according to method described
% in Benjamini & Hochberg, 1995, J Royal Stat Soci B, 57(1):289-300 and
% used in Guo et al (2014) Neuron
% Briefly, it ranks p-values in ascending order and defines a value i as being
% significant if it satisfies P(i) <= alpha*i/n where n is the number of
% comparisons
% alpha_correct is FDR-corrected statistical threshold
% isSig is nlocations x 1 boolean vector where true indicates significance
%
% LP aug 2016

if nargin < 2
  alpha = 0.05;
end

n             = numel(pvals);
[pranked,ix]  = sort(pvals,'ascend');

fdr           = ((1:n)*alpha)./n;
isSig         = pranked<=fdr';

if n == 1
    alpha_correct = alpha;
else
    if sum(isSig) > 0
        alpha_correct = pranked(sum(isSig));
    else
        alpha_correct = alpha;
    end
end

[~,nix]       = sort(ix,'ascend');
isSig         = isSig(nix);