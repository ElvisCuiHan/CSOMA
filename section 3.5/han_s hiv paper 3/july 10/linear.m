function result = linear(x, k)
    % 3 factors choose 2
    design = x(1:4*k);
    factor_choose = x(4*k+1:4*k+3) > 0.5;
    zero_index = find(~factor_choose);
    if abs(sum(factor_choose) - 2) ~= 0
        result = 100000000;
    else
        X = cat(2, ones(k, 1), reshape(design, [k, 4]));
        X(:, zero_index+1) = [];
        Info = zeros(3, 3);
        total_weight = sum(X(:, 4));
        X(:, 4) = X(:, 4) ./ total_weight;
        for i=1:k
            Info = Info + X(i, 4) * X(i, 1:3)' * X(i, 1:3);
        end
        result = -log(det(Info));
    
    end
    
end
