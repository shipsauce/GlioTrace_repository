function out = logsumexp_vec(v)
% v: column or row vector
m = max(v);
out = m + log(sum(exp(v - m)));
end
