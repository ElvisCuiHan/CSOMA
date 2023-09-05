%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Program:    Bayesian Optimal Design for Nonlinear Mixed Models Applied
%               to HIV Dynamics Using CSO-MA
%   Author:     Zizhao Zhang, Elvis Han Cui
%               hangzizhao.zzz@alibaba-inc.com
%               elvishancui@g.ucla.edu
%   Purpose:    Program calculates the Bayesian optimal design for a
%               nonlinear model applied to HIV dynamics suggested in Hans
%               and Chaloner (2004). The main optimization tool is
%               competitve swarm optimizer with mutated agents (CSO-MA).
%   Output:     design:     A 16x1 vector representing the design time
%                           points for the study.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% N: Number of iid samples used for each theta_i.
% K: Number of Gaussian quadratures.
% S: Number of different y-samples to approximate the expectation.
N = 5;
K = 5;
S = 5;

% Below are hyper-parameters for Wishart, Gaussian and Gamma distributions.
Omega = [0.26, 0, 0; 0, 2.5, 0; 0, 0, 2.5];
gamma_ = 3;
eta = [11, 1.1, -1];
Lambda = [6.0, 0, 0; 0, 0.1, 0; 0, 0, 0.01];
alpha = 4.5;
beta = 9.0;

% Below are parameters for the Markov chain.
num_mu = 4;
burnin = 20 ;
thinning = 2;

% Below are some candidate designs in Hans and Chaloner (2004) for
% comparison purposes.

%t1 = [0, 0, 0.083, 0.167, 0.417, 0.667, 0.917, 1.167, 1.417, 1.667, 1.917, 2.917, 3.917, 4.917, 5.917, 6.917];
%t2 = [0, 0, 0.083, 0.167, 0.417, 0.667, 0.917, 1.167, 1.417, 1.667, 1.917, 2, 2.917, 3, 6.834, 6.917];
%t5 = [0, 0, 0.021, 0.042, 0.083, 0.125, 0.167, 0.417, 0.667, 0.917, 1.167, 1.417, 1.667, 1.917, 2.917, 6.917];

% num_t: Length of the design (set to 16 in our paper).
% max_t: Maximum of designed time point (7 in our paper).
% upper_t: Upper bound for CSO-MA.
% lower_t: Lower bound for CSO-MA.
% swarmsize: Size of swarm (16 is sufficient to get reasonable results).
% phi: tuning parameter in CSO-MA (see https://link.springer.com/article/10.1007/s12293-020-00305-6).
% max_iter: Number of iterations for CSO-MA.

num_t = 16;
max_t = 7;
upper_t = max_t * ones(1, num_t);
lower_t = 0 * ones(1, num_t);
swarmsize = 16;
phi = 0;
max_iter = 30;

% obj: Objective function defined in our paper, i.e., the expectation of
%      posterior variance.
obj = @(t)F_3(t, S, N, K, Omega, gamma_, eta, Lambda, alpha, beta, num_mu, burnin, thinning);

% Let's do it!
tic
[fval, design] = csoma(obj, lower_t, upper_t, swarmsize, phi, max_iter);
toc
 
x0 = lower_t + unifrnd(1, num_t) .* (upper_t - lower_t);
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
fmincon(obj, x0, [], [], [], [], lower_t, upper_t, [], options);

