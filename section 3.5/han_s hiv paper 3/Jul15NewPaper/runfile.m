theta1 = [1.5 -0.7 0.4 -1.1 1.9 2.03 -1.14]';
k = 12;
interval = 5;
step = 0.5;
lb = cat(2, -interval*ones(1, 3*k), zeros(1, k));
ub = cat(2, interval*ones(1, 3*k), ones(1, k));
swarmsize = 400;
maxiter = 15000;
phi = 0;
%CLL = @(x)cll(x, theta1, k);
%[cf1, cd1] = cso(CLL, lb, ub, swarmsize, phi, maxiter);

%c11 = [];
%for x1 = -interval:step:interval
%    for x2 = -interval:step:interval
%        for x3 = -interval:step:interval
%            c11 = [c11 clleq(x1, x2, x3, cd1, theta1, k)];
%        end
%    end
%end
%max(c11)

%for i = 1:length(c11)
%    plot([i, i], [-8, c11(i)])
%    hold on
%end



LOGIS = @(x)logis(x, theta1, k);
[lf1, ld1] = cso(LOGIS, lb, ub, swarmsize, phi, maxiter);

l11 = [];
for x1 = -interval:step:interval
    for x2 = -interval:step:interval
        for x3 = -interval:step:interval
            l11 = [l11 logiseq(x1, x2, x3, ld1, theta1, k)];
        end
    end
end
max(l11)

reshape(ld1, [k, 4])

for i = 1:length(l11)
    plot([i, i], [-8, l11(i)])
    hold on
end


maxiter = 3000;
POI = @(x)poi(x, theta1, k);
[pf1, pd1] = cso(POI, lb, ub, swarmsize, phi, maxiter);

p11 = [];
for x1 = -interval:step:interval
    for x2 = -interval:step:interval
        for x3 = -interval:step:interval
            p11 = [p11 poieq(x1, x2, x3, pd1, theta1, k)];
        end
    end
end
max(p11)
reshape(pd1, [k, 4])
for i = 1:length(p11)
    plot([i, i], [-9, p11(i)])
    hold on
end
