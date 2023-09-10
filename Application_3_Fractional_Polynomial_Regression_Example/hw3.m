 
  %f2 = @(x)(det([0 0 1]*inv([1 sum(x(1:3).*x(4:6)/sum(x(4:6))) sum(x(1:3).^2.*x(4:6)/sum(x(4:6))); ...
  %    sum(x(1:3).*x(4:6)/sum(x(4:6))) sum(x(1:3).^2.*x(4:6)/sum(x(4:6))) sum(x(1:3).^3.*x(4:6)/sum(x(4:6))); ...
  %    sum(x(1:3).^2.*x(4:6)/sum(x(4:6))) sum(x(1:3).^3.*x(4:6)/sum(x(4:6))) sum(x(1:3).^4.*x(4:6)/sum(x(4:6)))]+0.00001*eye(3))*[0 0 1]'));
  f2 = @(x)(det([1 0 0; 0 0 1]*inv([1 sum(x(1:3).*x(4:6)/sum(x(4:6))) sum(x(1:3).^2.*x(4:6)/sum(x(4:6))); ...
      sum(x(1:3).*x(4:6)/sum(x(4:6))) sum(x(1:3).^2.*x(4:6)/sum(x(4:6))) sum(x(1:3).^3.*x(4:6)/sum(x(4:6))); ...
      sum(x(1:3).^2.*x(4:6)/sum(x(4:6))) sum(x(1:3).^3.*x(4:6)/sum(x(4:6))) sum(x(1:3).^4.*x(4:6)/sum(x(4:6)))]+0.00001*eye(3))*[1 0 0; 0 0 1]'));
  lb = [-1 -1 -1 0 0 0];
  ub = [1 1 1 1 1 1];
  swarmsize = 100;
  maxiter = 200;
  phi = 0;
  [r x] = cso(f2, lb, ub, swarmsize, phi, maxiter)
  
  
  y = [-1 0 1 0.25 0.5 0.25];
  M = [1 sum(y(1:3).*y(4:6)) sum(y(1:3).^2.*y(4:6)); ...
      sum(y(1:3).*y(4:6)) sum(y(1:3).^2.*y(4:6)) sum(y(1:3).^3.*y(4:6)); ...
      sum(y(1:3).^2.*y(4:6)) sum(y(1:3).^3.*y(4:6)) sum(y(1:3).^4.*y(4:6))];
  %A = [0 0 1];
  A = [1 0 0; 0 1 0];
  f = @(x)([1 x x^2]*inv(M)*A'*inv(A*inv(M)*A')*A*inv(M)*[1 x x^2]');
 
  fplot(f, [-1, 1])
  
  
  m1 = [3 -1/3 5/3; -1/3 5/3 -1/3; 5/3 -1/3 5/3];
  
  mf1 = @(x)(x^2 * [0 0 1]*inv(m1)*[0 0 1]');
  fplot(mf1, [-1, 1])
  
  m2 = [1 -1 1; -1 1 -1; 1 -1 1];
  mf2 = @(x)([0 0 x]*inv(m2 + 0.000001*eye(3))*[0 0 x]');
  fplot(mf2, [-1, 1])
  
  % 2 
  
cl = @(x)cll1(x, [0.1 0.3]', 2)
d = 6800
lb = cat(2, -d*ones(1, 2), zeros(1, 2));
ub = cat(2, d*ones(1, 2), ones(1, 2));
swarmsize = 200;
maxiter = 800;
phi = 0;
[r x] = cso(cl, lb, ub, swarmsize, phi, maxiter)
  
  
  cl = @(x)cll2(x, [0.1 0.3 -1]', 4)
  lb = [-1 -1 -1 -1 0 0 0 0];
  ub = [1 1 1 1 1 1 1 1];
  swarmsize = 100;
  maxiter = 400;
  phi = 0;
  [r x] = cso(cl, lb, ub, swarmsize, phi, maxiter)
  
  
  % FED
f2 = @(x)(det([1 0 0; 0 0 1]*inv([1 sum(x(1:3).*x(4:6)/sum(x(4:6))) sum(x(1:3).^2.*x(4:6)/sum(x(4:6))); ...
      sum(x(1:3).*x(4:6)/sum(x(4:6))) sum(x(1:3).^2.*x(4:6)/sum(x(4:6))) sum(x(1:3).^3.*x(4:6)/sum(x(4:6))); ...
      sum(x(1:3).^2.*x(4:6)/sum(x(4:6))) sum(x(1:3).^3.*x(4:6)/sum(x(4:6))) sum(x(1:3).^4.*x(4:6)/sum(x(4:6)))]+0.00001*eye(3))*[1 0 0; 0 0 1]'));
  
options = optimoptions('particleswarm','MinNeighborsFraction',1);
options.SwarmSize = 20;
%options.inertia = 1;
options.SelfAdjustmentWeight = 2;
options.SocialAdjustmentWeight = 2;
options.MaxIterations = 20;
lb = [-1 -1 -1 0 0 0];
ub = [1 1 1 1 1 1];
nvars = 6;
tic
[x, fval, exitflag] = particleswarm(f2, nvars, lb, ub, options)
toc
  

% 1b PSO 20 iterations 400 function evaluations 0.01s

cl = @(x)cll2(x, [0.1 0.3 -1]', 4)
%lb = [-1000 -100 -100 -100 0 0 0 0];
lb = cat(2, -10000*ones(1, 4), zeros(1, 4));
ub = cat(2, 10000*ones(1, 4), ones(1, 4));
options = optimoptions('particleswarm','MinNeighborsFraction',1);
options.SwarmSize = 120;
%options.inertia = 1;
options.SelfAdjustmentWeight = 2;
options.SocialAdjustmentWeight = 2;
options.MaxIterations = 300;
nvars = 8;
tic
[x, fval, exitflag] = particleswarm(cl, nvars, lb, ub, options)
toc

% 2b pso 20 iterations 400 function evaluations 0.02s

% 3b pso 30 iterations 900 function evaluations 0.08s

bd = [-1 -1 -1 0 0 0] + rand(1, 6) .* [2 2 2 1 1 1];
options = optimoptions('particleswarm','MinNeighborsFraction',1);
options.SwarmSize = 10;
%options.inertia = 1;
options.SelfAdjustmentWeight = 2;
options.SocialAdjustmentWeight = 2;
options.MaxIterations = 20;
options.OptimalityTolerance = 1e-4;
tic
fval = 100000;
numoffe = 0;
for iter = 1:4
    for j = 1:3
        if j == 1
            f = @(x)f2([x(1) bd(2) bd(3) x(2) bd(5) bd(6)]);
        elseif j == 2
            f = @(x)f2([bd(1) x(1) bd(3) bd(4) x(2) bd(6)]);
        else
            f = @(x)f2([bd(1) bd(2) x(1) bd(4) bd(5) x(2)]);
        end
        lb = [-1 0];
        ub = [1 1];
        nvar = 2;
        [x, fvalt, exitflag, output] = particleswarm(f, nvar, lb, ub, options)
        if j == 1
            bd = [x(1) bd(2) bd(3) x(2) bd(5) bd(6)];
        elseif j == k
            bd = [bd(1) x(1) bd(3) bd(4) x(2) bd(6)];
        else
            bd = [bd(1) bd(2) x(1) bd(4) bd(5) x(2)];
        end
        if fval - fvalt < 0.0001
            break;
        end
    end 
end
toc
bd
  




% 1b PSO 20 iterations 400 function evaluations 0.01s

qq = @(x)q5(x, [4.298 0.05884 21.8]', 4);
lb = [0 0 0 0 0 0 0 0];
ub = [100 100 100 100 1 1 1 1];
options = optimoptions('particleswarm','MinNeighborsFraction',1);
options.SwarmSize = 40;
%options.inertia = 1;
options.SelfAdjustmentWeight = 2;
options.SocialAdjustmentWeight = 2;
options.MaxIterations = 200;
nvars = 8;
tic
[x, fval, exitflag] = particleswarm(qq, nvars, lb, ub, options)
toc
  
n = 100;
theta = zeros(3, n);
theta(1, :) = unifrnd(0.04884, 0.06884, [1 n]);
theta(2, :) = unifrnd(3.298, 5.298, [1 n]);
theta(3, :) = 21.8 * ones(1, n);
qq = @(x)q5(x, theta, 4);
lb = [0 0 0 0 0 0 0 0];
ub = [100 100 100 100 1 1 1 1];
options = optimoptions('particleswarm','MinNeighborsFraction',1);
options.SwarmSize = 40;
%options.inertia = 1;
options.SelfAdjustmentWeight = 2;
options.SocialAdjustmentWeight = 2;
options.MaxIterations = 200;
nvars = 8;
tic
[x, fval, exitflag] = particleswarm(qq, nvars, lb, ub, options)
toc
  