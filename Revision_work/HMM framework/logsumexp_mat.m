function out = logsumexp_mat(M, dim)
% M: 2D matrix
% dim = 1  -> reduce over rows, output is 1 x N
% dim = 2  -> reduce over cols, output is M x 1

m = max(M, [], dim);
M_shift = M - m;
out = m + log(sum(exp(M_shift), dim));
end
