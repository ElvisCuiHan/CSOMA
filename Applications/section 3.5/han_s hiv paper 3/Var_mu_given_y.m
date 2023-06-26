function var_mu_1_given_y = Var_mu_given_y(t, Omega, gamma_, eta, Lambda, alpha, beta, l, m)
   %{
                                                                           -> (t_1) y_{111} 
                                                                           -> (t_2) y_{112}    
                                   -> theta_{111}, theta_{112}, theta_{113} ...
                                                                           -> (t_n) y_{11n}
                                   ...
         mu_{11}, mu_{12}, mu_{13}  
                                                                           -> (t_1) y_{1m1} 
                                                                           -> (t_2) y_{1m2}    
                                   -> theta_{1m1}, theta_{1m2}, theta_{1m3} ...
                                                                           -> (t_n) y_{1mn}
        .
        .
        .
        
                                                                           -> (t_1) y_{111} 
                                                                           -> (t_2) y_{l12}    
                                   -> theta_{l11}, theta_{l12}, theta_{l13} ...
                                                                           -> (t_n) y_{l1n}
                                   ...
         mu_{l1}, mu_{l2}, mu_{l3}  
                                                                           -> (t_1) y_{lm1} 
                                                                           -> (t_2) y_{lm2}    
                                   -> theta_{lm1}, theta_{lm2}, theta_{lm3} ...
                                                                           -> (t_n) y_{lmn}


   %}
   
   n = length(t);
   p = 3;
   
   
   mu = mvnrnd(eta, Lambda, l);
   
   %mu = zeros(k*n, p);
   theta = zeros(l*m, p);
   s_theta = zeros(1*m*n, 1);
   y = zeros(l*m*n, 1);
   
   for i = 1:l
       Sigma_inv = wishrnd(Omega, gamma_);
       sigma_inv_sq = gamrnd(alpha, beta);
       %mu(i*n+1:(i+1)*n, :) = repmat(mvnrnd(eta, Lambda), n, 1);
       theta((i-1)*m+1:i*m, :) = mvnrnd(mu(i, :), inv(Sigma_inv), m);
       log_V = theta((i-1)*m+1:i*m, 1);
       
       c = exp(theta((i-1)*m+1:i*m, 2));
       
       delta = exp(theta((i-1)*m+1:i*m, 3));
       
       s = log_V + log(c.^2.*exp(-delta.*t)./(c-delta).^2 - (c.^2-(c-delta).^2).*exp(-c.*t)./(c-delta).^2 - c.*delta.*t.*exp(-c.*t)./(c-delta));
       size(s)
       s_theta((l-1)*m*n+1:l*m*n) = reshape(s', 1, n*m);
       y((l-1)*m*n+1:l*m*n) = normrnd(s_theta((l-1)*m*n+1:l*m*n), sqrt(1/sigma_inv_sq));
   end
   
   
   constant = beta^alpha * gamma(alpha+0.5) / (2 * pi)^2 / gamma(alpha) / 2^gamma_ / pi^1.5 / gamma((gamma_-1)/2) ...
            / gamma(gamma_/2-1) / (det(Omega))^(gamma_/2);
        
   f_mu = 1 / (2*pi)^(3/2) / det(Lambda) .* exp(-sum((mu - eta) * inv(Lambda) .* (mu - eta), 2)); % 1 * l
 
   
   Omega_inv = inv(Omega);
   
   integration_theta_long_1 = 1 ./ ((y - s_theta).^2 + 2 * beta).^(alpha + 0.5);
   
   integration_theta_long_2 = [];
   
   
   for i = 1:l*m
       tmp_matrix = Omega_inv + (theta(i, :) - mu(ceil(i/m), :))' * (theta(i, :) - mu(ceil(i/m), :));
       integration_theta_long_2 = [integration_theta_long_2 (det(inv(tmp_matrix))).^(gamma_/2)];
   end
   integration_theta_long = integration_theta_long_1 + repelem(integration_theta_long_2', n); % 1 * lm
   integration_theta = accumarray(ceil((1:numel(integration_theta_long))/(m*n))',integration_theta_long(:),[],@mean); % 1 * l
   
 
   f_mu_given_y = constant * f_mu .* integration_theta / length(integration_theta);
   
   mu_1_sq_mean = (mu(:, 1).^2)' * f_mu_given_y;
   
   mu_1_mean = mu(:, 1)' * f_mu_given_y / length(f_mu_given_y);
   
   var_mu_1_given_y = mu_1_sq_mean - mu_1_mean^2 ;
end
 
