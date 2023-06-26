function result = cll(x, theta, k)
    digits = 2000;
    num_theta = size(theta, 2);
    PX = cat(2, ones(k, 1), reshape(x, [k, 3]));
    X = [PX(:, 1:3) PX(:, 2).*PX(:, 3) PX(:, 4)];
    long_Info = zeros(4, num_theta);
    total_weight = sum(X(:, 5));
    X(:, 5) = X(:, 5) ./ total_weight;
    linear_part = X(:, [1:4]) * theta;
    exp_part = exp(-exp(linear_part)) .* exp(linear_part) .* exp(linear_part) ./ (1 - exp(-exp(linear_part)));
    for i=1:4
        for j=1:4
            long_Info((i-1)*4+j, :) = X(:, 5)' .* X(:, i)' .* X(:, j)' * exp_part;
        end
    end
    matrix_Info = reshape(long_Info, [4, 4, num_theta]);
    result = mean(arrayfun(@(K) -log(det(matrix_Info(:, :, K))), 1:num_theta));
end
