function result = frac1_eq(x, design, k, D, m)
    
    % 1, x_i^(-1/3), x_i^(-0.5), x_i^(1/3), x_i^0.5, x_i, x_i^2
    xpart = design(1:k*m);
    wpart = design(k*m+1:k*m+k);
    PX = cat(2, ones(k*m, 1), reshape(xpart, [k*m, 1])); %only xs and weights
    X = [PX(:, 1) PX(:, 2).^(-0.333) PX(:, 2).^(-0.5) PX(:, 2).^(0.333) PX(:, 2).^(0.5) PX(:, 2) PX(:, 2).^2];
    M = zeros(7, 7);
    total_weight = sum(wpart);
    wpart = wpart ./ total_weight;
    for i = 0:(k-1)
        tmp = X(i*m+1:(i+1)*m, 1:7) ;
        M = M + tmp' / (eye(m) + tmp * D * tmp') * tmp * wpart(i+1);
    end
    
    tt = [ones(m, 1), x'.^(-0.333), x'.^(-0.5), x'.^(0.333), x'.^0.5, x', x'.^2];
    
  
    result = 7 - trace((tt' / (eye(m) + tt * D * tt') * tt) / M);
    
    
end
