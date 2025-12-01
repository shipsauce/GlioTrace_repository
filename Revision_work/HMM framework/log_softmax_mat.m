function out = log_softmax_mat(x, dim)
% Log softmax along dimension dim
m = max(x, [], dim);
x_shift = x - m;
out = x_shift - log(sum(exp(x_shift), dim));
end
