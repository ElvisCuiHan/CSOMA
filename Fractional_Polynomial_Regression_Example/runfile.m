k = 4;
region_ub = [-1, 2, 5];
region_lb = [-5, -1, 1];
lb = cat(2, region_lb(1)*ones(1, k), region_lb(2)*ones(1, k), region_lb(3)*ones(1, k), zeros(1, k+3));
ub = cat(2, region_ub(1)*ones(1, k), region_ub(2)*ones(1, k), region_ub(3)*ones(1, k), ones(1, k+3));
swarmsize = 400;
maxiter = 150 ;
phi = 0;

M = @(x)linear(x, k);

swarmsize = 64; % swarmsize
phi = 0; % parameter phi 
maxiter = 200; % iteration number 
tic
[fitness1, design1] = cso(M, lb, ub, swarmsize, phi, maxiter); % to find the optimal design
toc

design1(3*k+1:4*k) = design1(3*k+1:4*k) / sum(design1(3*k+1:4*k));
design1(4*k+1:4*k+3) = design1(4*k+1:4*k+3) > 0.5;
idx = find(~design1(4*k+1:4*k+3));
a = reshape(design1(1:4*k), [k, 4]);
a(:, idx) = [];
a
fprintf("choose columns %d, %d \n", find(design1(4*k+1:4*k+3)))

 

%%%%%%%%%%%%%%%%%%%%%%%%%%
D = diag([0.1, 0.2, 0.4, 0.5, 0.5, 0.3, 0.2])
D1 = [1, 0.8, 0.4; 0.8, 1.2, 0.5; 0.4, 0.5, 0.9];
D2 = [1.3, 0.6, 0.6, 0.3; 0.6, 1.2, 0.7, 0.4; 0.6, 0.7, 1.3, 0.3; 0.3, 0.4, 0.3, 1];
D = blkdiag(D1, D2)
m = 4;
k = 8; 
lb = cat(2,  ones(1, k), zeros(1, k));
ub = cat(2,  10 * ones(1, k), ones(1, k));
swarmsize = 3200;
maxiter = 2500 ;
phi = 0.05;

M = @(x)frac(x, k, D, m);

 
tic
[fitness1, design1] = cso(M, lb, ub, swarmsize, phi, maxiter); % to find the optimal design
toc

design1(k+1:2*k) = design1(k+1:2*k) / sum(design1(k+1:2*k));
reshape(design1, [k, 2])



eq = @(x)frac_eq(x, design1, k, D, m);
lb = 1;
ub = 10;
swarmsize = 100;
maxiter = 150;
phi = 0;
tic
[f, fx] = cso(eq, lb, ub, swarmsize, phi, maxiter); % to find the optimal design
toc


eqp = @(x)(0 - eq(x))
fplot(eqp, [lb, ub])



%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
D = diag([0.1, 0.2, 0.4, 0.5, 0.5, 0.3, 0.2])
D1 = [1, 0.8, 0.4; 0.8, 1.2, 0.5; 0.4, 0.5, 0.9];
D2 = [1.3, 0.6, 0.6, 0.3; 0.6, 1.2, 0.7, 0.4; 0.6, 0.7, 1.3, 0.3; 0.3, 0.4, 0.3, 1];
D = blkdiag(D1, D2)
m = 4;
k = 8; 
lb = cat(2,  ones(1, k*m), zeros(1, k));
ub = cat(2,  10 * ones(1, k*m), ones(1, k));
swarmsize = 3200;
maxiter = 1500 ;
phi = 0.05;

M = @(x)frac1(x, k, D, m);

 
tic
[fitness1, design1] = cso(M, lb, ub, swarmsize, phi, maxiter); % to find the optimal design
toc
 

w = design1(k*m+1:k*m+k) / sum(design1(k*m+1:k*m+k));
d = reshape(design1(1:k*m), [m, k])'
cat(2, d, w')
 


eq = @(x)frac1_eq(x, design1, k, D, m);
lb = ones(1, m);
ub = 10 * ones(1, m);
swarmsize = 200;
maxiter = 350;
phi = 0;
tic
[f, fx] = cso(eq, lb, ub, swarmsize, phi, maxiter); % to find the optimal design
toc


ee = [];
for t1 = 1:1:10
    for t2 = 1:1:10    
        for t3 = 1:1:10
            for t4 = 1:1:10
                tmp = 0 - eq([t1, t2, t3, t4]);
                 
                ee = [ee, tmp-0.1];
            end
        end
    end
end

plot(ee)

eqp = @(x)(0 - eq(x))
fplot(eqp, [lb, ub])




