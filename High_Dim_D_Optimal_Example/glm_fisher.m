function result = glm_fisher(b, theta)
    % b:     6Mx1 parameter
    % theta: 16x1 vector of nominal values
    b = b';
    %disp(size(b));
    M = size(b, 1) / 6;
    cache = zeros(size(theta, 1), size(theta, 1));
    %disp(M);
    for i = 1:M
        p(i) = b(6*i);
    end
    p = (p + 1);
    p = p / sum(p);

    for i = 1:M
        bi = b(6*i-5:6*i-1);
        pi = p(i);
        %disp(pi);
        cache = cache + single_fisher(bi, pi, theta);
    end
    %disp(cache);
    result = - logdet(cache);
    %result = -log(det(cache));
end

function sf = single_fisher(b, p, theta)
    % b:     5x1 parameter
    % p:     1x1 parameter
    % theta: 16x1 vector of nominal values
    x1 = b(1); 
    x2 = b(2);
    x3 = b(3);
    x4 = b(4);
    x5 = b(5);
    x = [1 x1 x2 x3 x4 x5 x1*x2 x1*x3 x1*x4 x1*x5 x2*x3 x2*x4 x2*x5 ...
        x3*x4 x3*x5 x4*x5];
    eta = x * theta;
    %w = logit_weight(eta);
    w = poisson_reg_weight(eta);
    sf = (p * w) * (x' * x);
end

function omega = logit_weight(eta)
    % eta: canonical link
    omega = exp(eta) / (1 + exp(eta)) / (1 + exp(eta));
end

function omega = poisson_reg_weight(eta)
    % eta: canonical link
    omega = exp(eta);
end

function y = logdet(A)
% log(det(A)) where A is positive-definite.
% This is faster and more stable than using log(det(A)).

% Written by Tom Minka
% (c) Microsoft Corporation. All rights reserved.

U = chol(A);
y = 2*sum(log(diag(U)));
end