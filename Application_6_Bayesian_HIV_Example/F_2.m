function res = F_2(t, y, mu, N, K, Omega, gamma_, eta, Lambda, alpha, beta)
   % y: 1 * T
   % t: 1 * T
   % mu: 1 * T
   
   tmp_res = 0;
   
   Omega_chol = chol(Omega);
   
   f_mu = mvnpdf(mu, eta, Lambda);
   
   for i = 1:N
       Sigma_inv = wishrnd(Omega, gamma_, Omega_chol);
       sigma_inv_sq = gamrnd(alpha, beta);
      
       Sigma = inv(Sigma_inv + 0.000001 * eye(3));
       sigma_sq = 1 / sigma_inv_sq;
       
       F_1_res = F_1(t, y, mu, Sigma, sigma_sq, K);
  
       tmp_res = tmp_res + F_1_res;
   end
   
   res = f_mu * tmp_res / N;
end 
 
