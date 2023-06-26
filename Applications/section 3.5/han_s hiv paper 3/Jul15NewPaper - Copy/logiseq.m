function result = logiseq(x1, x2, design, theta, k)
    digits = 1200;
    num_theta = size(theta, 2);
    PX = cat(2, ones(k, 1), reshape(design, [k, 3]));
    X = [PX(:, 1:3) PX(:, 2).*PX(:, 3) PX(:, 4)];
    long_Info = zeros(4, num_theta);
    total_weight = sum(X(:, 5));
    X(:, 5) = X(:, 5) ./ total_weight;
    linear_part = X(:, [1:4]) * theta;
    exp_part = exp(linear_part) ./ (1 + exp(linear_part)).^2;
    for i=1:4
        for j=1:4
            long_Info((i-1)*4+j, :) = X(:, 5)' .* X(:, i)' .* X(:, j)' * exp_part;
        end
    end
    matrix_Info = reshape(long_Info, [4, 4, num_theta]);
    Info_inv = arrayfun(@(K) [1 x1 x2 x1*x2] / matrix_Info(:, :, K) ...
               * [1 x1 x2 x1*x2]', 1:num_theta);
    result = mean(exp([1 x1 x2 x1*x2] * theta)./ ...
               (1 + exp([1 x1 x2 x1*x2] * theta)).^2 .* Info_inv) - 4;
end
