Omega = [0.26, 0, 0; 0, 2.5, 0; 0, 0, 2.5];
gamma_ = 3;
eta = [11, 1.1, -1];
Lambda = [6.0, 0, 0; 0, 0.1, 0; 0, 0, 0.01];
alpha = 4.5;
beta = 9.0;



%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
% Start with a sequence of random y;
% Use cso-ma to find best tl;
% Given that best t, we resample y;
% Iterate step 2 and 3;

k = 100;  % k: size for theta
ss = 200; % ss: size for mu, Sigma, sigma_sq
max_epoch = 50; 
num_t = 2;
max_t = 3;
upper_t = max_t * ones(1, num_t);
lower_t = 0 * ones(1, num_t);
swarmsize = 20;
phi = 0;
max_iter = 20;

sample_hyper = 100;
num_y = sample_hyper * num_t;
% y = unifrnd(-5, 5, num_y);

Tm = zeros(max_epoch, num_t);

fres = [];

t = [0.5, max_t];
for epoch = 1:max_epoch
    
    % After getting the temporary optimal design, we resample y.
    y = [];
    for j = 1:sample_hyper
        Sigma_inv = wishrnd(Omega, gamma_); % 3 * 3
        sigma_inv_sq = gamrnd(alpha, beta); % a scalar
        mu = mvnrnd(eta, Lambda); % 1 * 3
        theta = mvnrnd(mu, inv(Sigma_inv)); % 1 * 3
    
  
        log_V = theta(1);  % a scalar
       
        c = exp(theta(2));  % a scalar
        
        delta = exp(theta(3));  % a scalar
       
        s = log_V + log(c.^2.*exp(-delta.*t)./(c-delta).^2 - (c.^2-(c-delta).^2).*exp(-c.*t)./(c-delta).^2  ...
            - c.*delta.*t.*exp(-c.*t)./(c-delta));  % 1 * length(t);
        
        tmp_y = normrnd(s, 1/sigma_inv_sq, 1, num_t);
        tmp_y = tmp_y .* (10000 >= tmp_y) + 10000 .* (1000 < tmp_y);
        tmp_y = tmp_y .* (-1000 < tmp_y) - 1000 .* (-1000 >= tmp_y);
        y = [y tmp_y];
    end
    
    obj = @(t)F_3(t, y, k, Omega, gamma_, eta, Lambda, alpha, beta, ss);

    tic
    [fval, design] = csoma(obj, lower_t, upper_t, swarmsize, phi, max_iter);
    toc
    
    t = design;
    fres = [fres fval];
    fval
    Tm(epoch, :) = sort(t);
end

for i = 1:num_t
    plot(Tm(:, i), 'o-')
    hold on
end


fres

