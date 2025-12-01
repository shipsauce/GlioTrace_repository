function A = calc_A_global(xis)

K = numel(xis);
N = size(xis{1}, 2);

A_num = zeros(N, N);
A_den = zeros(N, 1);

for k = 1:K
    xi = xis{k};            % [N x N x (T-1)]

    A_num = A_num + squeeze(sum(xi, 1));                  % sum over time
    A_den = A_den + squeeze(sum(sum(xi, 3), 1)).';        % sum_j xi_t(i,j)
end

A = A_num ./ A_den;  % row normalize

end