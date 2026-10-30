function selected = weighted_sample_without_replacement(weights, n, algorithm)
%WEIGHTED_SAMPLE_WITHOUT_REPLACEMENT Select distinct indices by weight.
% 'legacy' preserves the published implementation's floating-point behavior.
% 'log_keys' avoids the numerical underflow of the original power keys.
if nargin < 3 || isempty(algorithm)
    algorithm = 'legacy';
end
weights = weights(:);
nTotal = numel(weights);
n = min(max(0,round(n)),nTotal);
if n == 0
    selected = zeros(0,1);
    return;
end
weights(isnan(weights)) = 0;
finiteWeights = weights(isfinite(weights));
if isempty(finiteWeights)
    selected = randperm(nTotal,n)';
    return;
end
weights(isinf(weights)) = max(finiteWeights);
weights(weights < 0) = 0;
if sum(weights) <= 0
    selected = randperm(nTotal,n)';
    return;
end
weights = weights ./ sum(weights);
% Draw exactly one uniform random variate per candidate in both modes.
u = rand(nTotal,1);
positive = weights > 0;
if strcmpi(char(algorithm),'legacy')
    keys = zeros(nTotal,1);
    keys(positive) = u(positive).^(1 ./ weights(positive));
    keys(~positive) = -inf;
elseif strcmpi(char(algorithm),'log_keys')
    keys = -inf(nTotal,1);
    keys(positive) = log(max(u(positive),realmin)) ./ weights(positive);
else
    error('VGPP:SamplerMode','Unknown weighted sampler: %s',char(algorithm));
end
[~,order] = sort(keys,'descend');
selected = order(1:n);
end
