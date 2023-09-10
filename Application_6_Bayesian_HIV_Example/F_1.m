function res = F_1(t, y, mu, Sigma, sigma_sq, K)
    % y: 1 * T
    % t: 1 * T
    % mu: 1 * 3
    % Sigma_inv: 3 * 3

    Sigma = (Sigma + transpose(Sigma)) / 2;
    theta = mvnrnd(mu, Sigma, K^3);   % K * 3
    %theta = md_gauss(K, 3) * chol(Sigma) + repmat(mu, K^3, 1);

    
    T = length(t);

    log_V = theta(:, 1);  % K * 1
       
    c = exp(theta(:, 2));  % K * 1
       
    delta = exp(theta(:, 3));  % K * 1
       
    s = log_V + log(c.^2.*exp(-delta.*t)./(c-delta).^2 - (c.^2-(c-delta).^2).*exp(-c.*t)./(c-delta).^2  ...
        - c.*delta.*t.*exp(-c.*t)./(c-delta));  % K * T
    
    part = -sum((y - s).^2 / 2 / sigma_sq, 2); % K * 1
    mu_copy = repmat(mu, K^3, 1); % K * 3
    part2 = -0.5 .* diag((theta - mu_copy) * inv(Sigma) * (theta - mu_copy)'); % K * 1
    part = part + part2;
    
    const = 1 / (2 * pi)^(T/2) / sqrt(sigma_sq)^T / K;
    res = sum(exp(part)) * const;
    
end
 
