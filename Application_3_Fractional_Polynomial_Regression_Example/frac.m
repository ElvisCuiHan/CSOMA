function result = frac(x, k, D, m)
    
    % 1, x_i^(-1/3), x_i^(-0.5), x_i^(1/3), x_i^0.5, x_i, x_i^2
    PX = cat(2, ones(k, 1), reshape(x, [k, 2])); %only xs and weights
    X = [PX(:, 1) PX(:, 2).^(-0.333) PX(:, 2).^(-0.5) PX(:, 2).^(0.333) PX(:, 2).^(0.5) PX(:, 2) PX(:, 2).^2 PX(:, 3)];
    M = zeros(7, 7);
    total_weight = sum(X(:, 8));
    X(:, 8) = X(:, 8) ./ total_weight;
    for i = 1:k
        tmp = X(i, 1:7) .* ones(m, 7);
        M = M + tmp' / (eye(m) + tmp * D * tmp') * tmp * X(i, 8);
    end
    
    
    result = -log(det(M));
    
    
end
