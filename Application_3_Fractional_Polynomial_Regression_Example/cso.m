
function [minf, minx] = cso(obj_fun, lb, ub, swarmsize, phi, maxiter)
    assert(length(lb) == length(ub), 'Not equal length of bounds');
    if all(ub - lb <= 0) > 0
        error('Error. \n Upper bound must be greater than lower bound.')
    end
    vhigh = abs(ub - lb);
    vlow = -vhigh;
    S = swarmsize;
    D = length(ub);
    x = rand(S, D);
    x = bsxfun(@plus, lb, bsxfun(@times, ub-lb, x)); % randomly initalize all particles
    v = zeros([S D]); % set initial velocities to 0
    iter = 0;
    pairnum_1 = floor(S / 2);
    losers = 1:S;
    fx = arrayfun(@(K) obj_fun(x(K, :)), 1:S);
    randperm_index = randperm(S);
    while iter <= maxiter
        fx(losers) = arrayfun(@(K) obj_fun(x(K, :)), losers);
        swarm_center = mean(x); % calculate center all particles 
        randperm_index = randperm(S); % randomly permuate all particle indexes
        rpairs = [randperm_index(1:pairnum_1); randperm_index(S-pairnum_1+1:S)]'; % random pair
        cmask= (fx(rpairs(:, 1)) > fx(rpairs(:, 2)))'; 
        losers = bsxfun(@times, cmask, rpairs(:, 1)) + bsxfun(@times, ~cmask, rpairs(:, 2)); % losers who with larger values
        winners = bsxfun(@times, ~cmask, rpairs(:, 1)) + bsxfun(@times, cmask, rpairs(:, 2)); % winners who with smaller values
        R1 = rand(pairnum_1, D);
        R2 = rand(pairnum_1, D);
        R3 = rand(pairnum_1, D);
        v(losers, :) = bsxfun(@times, R1, v(losers, :)) + bsxfun(@times, R2, x(winners, :) - x(losers, :)) + phi * bsxfun(@times, R3, bsxfun(@minus, swarm_center, x(losers, :)));
        x(losers, :) = x(losers, :) + v(losers, :);
        maskl = bsxfun(@lt, x(losers, :), lb);
        masku = bsxfun(@gt, x(losers, :), ub);
        mask = bsxfun(@lt, x(losers, :), lb) | bsxfun(@gt, x(losers, :), ub);
        x(losers, :) = bsxfun(@times, ~mask, x(losers, :)) + bsxfun(@times, lb, maskl) + bsxfun(@times, ub, masku);
        iter = iter + 1;    
        fprintf('Iter: %d\n', iter);
        fprintf('Best fitness: %e\n', min(fx));
    end
    fprintf('Best fitness: %e\n', min(fx));
    [minf, min_index] = min(fx);
    minx = x(min_index, :);
    end
 
