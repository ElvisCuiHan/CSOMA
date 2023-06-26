function result = frac1(x, k, D, m)
    
    % 1, x_i^(-1/3), x_i^(-0.5), x_i^(1/3), x_i^0.5, x_i, x_i^2
    xpart = x(1:k*m);
    wpart = x(k*m+1:k*m+k);
    PX = cat(2, ones(k*m, 1), reshape(xpart, [k*m, 1])); %only xs and weights
    X = [PX(:, 1) PX(:, 2).^(-0.333) PX(:, 2).^(-0.5) PX(:, 2).^(0.333) PX(:, 2).^(0.5) PX(:, 2) PX(:, 2).^2];
    M = zeros(7, 7);
    total_weight = sum(wpart);
    wpart = wpart ./ total_weight;
    for i = 0:(k-1)
        tmp = X(i*m+1:(i+1)*m, 1:7) ;
        M = M + tmp' / (eye(m) + tmp * D * tmp') * tmp * wpart(i+1);
    end
    
    
    result = -log(det(M));
    
    
end
