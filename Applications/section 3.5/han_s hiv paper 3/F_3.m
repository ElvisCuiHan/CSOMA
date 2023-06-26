function res = F_3(t, M, N, K, Omega, gamma_, eta, Lambda, alpha, beta, num_mu, burnin, thinning)
   % t: 1 * T
   res_mu = zeros(M, num_mu);

   T = length(t);

   Omega_chol = chol(Omega);
   Sigma_inv = wishrnd(Omega, gamma_, Omega_chol); % 3 * 3
   sigma_inv_sq = gamrnd(alpha, beta); % a scalar
   mu_c = mvnrnd(eta, Lambda); % 1 * 3
   
   Sigma = inv(Sigma_inv + 0.000001 * eye(3));
   Sigma = (Sigma + transpose(Sigma)) / 2;
   
   stan_dev = sqrt(1 / sigma_inv_sq); 
   for i = 1:M
        theta = mvnrnd(mu_c, Sigma); % 1 * 3
        log_V = theta(1);  % a scalar
       
        c = exp(theta(2));  % a scalar
        delta = exp(theta(3));  % a scalar
       
        s = log_V + log(c.^2.*exp(-delta.*t)./(c-delta).^2 - (c.^2-(c-delta).^2).*exp(-c.*t)./(c-delta).^2  ...
            - c.*delta.*t.*exp(-c.*t)./(c-delta));  % 1 * T;
        
        y = normrnd(s, stan_dev, 1, T);
        
        mu = mu_c;
        % generate mu by y
        
        % 0.005
        % 0.0005
        cov_m = 0.0005 * eye(3);
        
        for j = 1:burnin
            z = mvnrnd(mu, cov_m);
            r = unifrnd(0, 1);
            tau = min(F_2(t, y, z, N, K, Omega, gamma_, eta, Lambda, alpha, beta) / F_2(t, y, mu, N, K, Omega, gamma_, eta, Lambda, alpha, beta), 1);
            if r <= tau
                mu = z;
            end
        end
        
        count = 1;
        flag = 1;
        while count <= num_mu
            z = mvnrnd(mu, cov_m);
            r = unifrnd(0, 1);
            tau = min(F_2(t, y, z, N, K, Omega, gamma_, eta, Lambda, alpha, beta) / F_2(t, y, mu, N, K, Omega, gamma_, eta, Lambda, alpha, beta), 1);
            if r <= tau
                mu = z;
            end
            if rem(flag, thinning) == 0
                res_mu(i, count) = mu(3);
                count = count + 1;    
            end 
            flag = 1 + flag;
        end

   end

   res = mean(var(res_mu, 0, 2));
end 
 
