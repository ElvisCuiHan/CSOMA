function result = logis(x, theta, k)
    digits = 1500;
    num_theta = size(theta, 2);
    PX = cat(2, ones(k, 1), reshape(x, [k, 4]));
    X = [PX(:, 1:4) PX(:, 2).*PX(:, 3) PX(:, 2).*PX(:, 4) PX(:, 3).*PX(:, 4) PX(:, 5)];
    long_Info = zeros(7, num_theta);
    total_weight = sum(X(:, 8));
    X(:, 8) = X(:, 8) ./ total_weight;
    linear_part = X(:, 1:7) * theta;
    exp_part = exp(linear_part) ./ (1 + exp(linear_part)).^2;
    for i=1:7
        for j=1:7
            long_Info((i-1)*7+j, :) = X(:, 8)' .* X(:, i)' .* X(:, j)' * exp_part;
        end
    end
    matrix_Info = reshape(long_Info, [7, 7, num_theta]);
    result = mean(arrayfun(@(K) -log(det(matrix_Info(:, :, K))), 1:num_theta));
end
