function result = poieq(x1, x2, x3, design, theta, k)
    digits = 1500;
    num_theta = size(theta, 2);
    PX = cat(2, ones(k, 1), reshape(design, [k, 4]));
    X = [PX(:, 1:4) PX(:, 2).*PX(:, 3) PX(:, 2).*PX(:, 4) PX(:, 3).*PX(:, 4) PX(:, 5)];
    long_Info = zeros(7, num_theta);
    total_weight = sum(X(:, 8));
    X(:, 8) = X(:, 8) ./ total_weight;
    linear_part = X(:, 1:7) * theta;
    exp_part = exp(linear_part);
    for i=1:7
        for j=1:7
            long_Info((i-1)*7+j, :) = X(:, 8)' .* X(:, i)' .* X(:, j)' * exp_part;
        end
    end
    matrix_Info = reshape(long_Info, [7, 7, num_theta]);
    Info_inv = arrayfun(@(K) [1 x1 x2 x3 x1*x2 x1*x3 x2*x3] / matrix_Info(:, :, K) ...
               * [1 x1 x2 x3 x1*x2 x1*x3 x2*x3]', 1:num_theta);
    result = mean(exp([1 x1 x2 x3 x1*x2 x1*x3 x2*x3] * theta) .* Info_inv) - 7;
end
