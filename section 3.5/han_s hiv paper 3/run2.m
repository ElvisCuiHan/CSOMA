
M = 5;
N = 5;
K = 5;

Omega = [0.26, 0, 0; 0, 2.5, 0; 0, 0, 2.5];
gamma_ = 3;
eta = [11, 1.1, -1];
Lambda = [6.0, 0, 0; 0, 0.1, 0; 0, 0, 0.01];
alpha = 4.5;
beta = 9.0;

num_mu = 4 ;
burnin = 20 ;
thinning = 2;

%t1 = [0, 0, 0.083, 0.167, 0.417, 0.667, 0.917, 1.167, 1.417, 1.667, 1.917, 2.917, 3.917, 4.917, 5.917, 6.917];

%t2 = [0, 0, 0.083, 0.167, 0.417, 0.667, 0.917, 1.167, 1.417, 1.667, 1.917, 2, 2.917, 3, 6.834, 6.917];
%F_3(t2, M, N, K, Omega, gamma_, eta, Lambda, alpha, beta, num_mu, burnin, thinning)

num_t = 16;
max_t = 7;
upper_t = max_t * ones(1, num_t);
lower_t = 0 * ones(1, num_t);
swarmsize = 16;
phi = 0;
max_iter = 30;

obj = @(t)F_3(t, M, N, K, Omega, gamma_, eta, Lambda, alpha, beta, num_mu, burnin, thinning);

tic
[fval, design] = csoma(obj, lower_t, upper_t, swarmsize, phi, max_iter);
toc
 
x0 = lower_t + unifrnd(1, num_t) .* (upper_t - lower_t);
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
fmincon(obj, x0, [], [], [], [], lower_t, upper_t, [], options);

%{

Omega = [0.26, 0, 0; 0, 2.5, 0; 0, 0, 2.5];
gamma_ = 3;
eta = [11, 1.1, -1];
Lambda = [6.0, 0, 0; 0, 0.1, 0; 0, 0, 0.01];
alpha = 4.5;
beta = 9.0;



Omega = [0.0078, 0, 0; 0, 0.075, 0; 0, 0, 0.075];
gamma_ = 100;
eta = [11, 1.1, -1];
Lambda = [0.6, 0, 0; 0, 0.01, 0; 0, 0, 0.001];
alpha = 81;
beta = 0.5;




num_t = 3;
max_t = 7;
upper_t = max_t * ones(1, num_t);
lower_t = 0 * ones(1, num_t);
swarmsize = 16;
phi = 0;
max_iter = 40;

obj = @(t)F_3(t, M, N, K, Omega, gamma_, eta, Lambda, alpha, beta, num_mu, burnin, thinning);

tic
[fval, design] = csoma(obj, lower_t, upper_t, swarmsize, phi, max_iter);
toc
 


%}

