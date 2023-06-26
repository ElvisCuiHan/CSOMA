k = 8;
swarmsize = 200;
maxiter = 1000;
phi = 0.05;

interval = 1;
step = 0.005;
lb = cat(2, -interval*ones(1, 2*k), zeros(1, k));
ub = cat(2, interval*ones(1, 2*k), ones(1, k));
LOGIS = @(x)logis(x, theta1, k);
[lf1, ld1] = cso(LOGIS, lb, ub, swarmsize, phi, maxiter);





for i = 1:length(l11)
    plot([i, i], [-3, l11(i)], 'black');
    hold on;
end
xlim([0 length(l11)]);
ylim([-3 1]);
%set(gca,'xtick',[]);
saveas(gcf,'log1.png')






POI = @(x)poi(x, theta1, k);
[pf1, pd1] = cso(POI, lb, ub, swarmsize, phi, maxiter);


for i = 1:length(p11)
    plot([i, i], [-4, p11(i)], 'black');
    hold on;
end
xlim([0 length(p11)]);
ylim([-4 1]);
%set(gca,'xtick',[]);
saveas(gcf,'poi1.png')










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


theta2 = [-1.7 -1 2 -1]';


LOGIS = @(x)logis(x, theta2, k);
[lf2, ld2] = cso(LOGIS, lb, ub, swarmsize, phi, maxiter);


for i = 1:length(l12)
    plot([i, i], [-3, l12(i)], 'black');
    hold on;
end
xlim([0 length(l12)]);
ylim([-3 1]);
%set(gca,'xtick',[]);
saveas(gcf,'log2.png')


POI = @(x)poi(x, theta2, k);
[pf2, pd2] = cso(POI, lb, ub, swarmsize, phi, maxiter);

for i = 1:length(p12)
    plot([i, i], [-4, p12(i)], 'black');
    hold on;
end
xlim([0 length(p12)]);
ylim([-4 1]);
%set(gca,'xtick',[]);
saveas(gcf,'poi2.png')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





theta3 = [-3 -2 3 1]';

LOGIS = @(x)logis(x, theta3, k);
[lf3, ld3] = cso(LOGIS, lb, ub, swarmsize, phi, maxiter);
% 1.8912


for i = 1:length(l13)
    plot([i, i], [-4, l13(i)], 'black');
    hold on;
end
xlim([0 length(l13)]);
ylim([-4 1]);
%set(gca,'xtick',[]);
saveas(gcf,'log3.png')

POI = @(x)poi(x, theta3, k);
[pf3, pd3] = cso(POI, lb, ub, swarmsize, phi, maxiter);


for i = 1:length(p13)
    plot([i, i], [-4, p13(i)], 'black');
    hold on;
end
xlim([0 length(p13)]);
ylim([-4 1]);
%set(gca,'xtick',[]);
saveas(gcf,'poi3.png')

l11 = [];
for x1 = -interval:step:interval
    for x2 = -interval:step:interval
        l11 = [l11 logiseq(x1, x2, ld1, theta1, k)]; 
    end
end
max(l11)


l12 = [];
for x1 = -interval:step:interval
    for x2 = -interval:step:interval
        l12 = [l12 logiseq(x1, x2, ld2, theta2, k)];
    end
end
max(l12)


l13 = [];
for x1 = -interval:step:interval
    for x2 = -interval:step:interval
       l13 = [l13 logiseq(x1, x2, ld3, theta3, k)];
    end
end
max(l13)

p11 = [];
for x1 = -interval:step:interval
    for x2 = -interval:step:interval
       p11 = [p11 poieq(x1, x2, pd1, theta1, k)];
    end
end
max(p11)

p12 = [];
for x1 = -interval:step:interval
    for x2 = -interval:step:interval
       p12 = [p12 poieq(x1, x2, pd2, theta2, k)];
    end
end
max(p12)

p13 = [];
for x1 = -interval:step:interval
    for x2 = -interval:step:interval
       p13 = [p13 poieq(x1, x2, pd3, theta3, k)];
    end
end
max(p13)



pd3(2*k+1:3*k) = pd3(2*k+1:3*k) / sum(pd3(2*k+1:3*k));
reshape(pd3, [k, 3])

Dl1 = [1 -0.3307 0.1416; ...
       0.1077 -1 0.1416; ...
       1 0.3307  0.2333; ...
       -0.7743 -1 0.2333; ...
       -1 1 0.25];
   
Dl2 = [-1 -0.246 0.2466; ...
        1 -1 0.25; ...
        -0.569 1 0.1284; ...
        0.869 1 0.2466; ...
        -1 0.7127 0.1284];
    
Dl3 = [0.3661 0.317 0.25; ...
       -1 -0.398 0.25;
       -1 1 0.25; ...
       1 1 0.25];
   
Dp1 = [1 -0.5 0.25; ... 
       0.3333 -1 0.25; ...
       1 -1 0.25; ...
       -1 1 0.25];
   
Dp2 = [0 1 0.25; ...
       -1 0.3333 0.2 5; ... 
       -1 1 0.25; ...
       1 -1 0.25];

Dp3 = [0.2361 0.382 0.25; ...
       -1 0 0.25; ...
       1 1 0.25; ...
       -1 1 0.25];
lf1 = 9.9618
lf2 = 10.9204
lf3 = 11.7832
pf1 = -1.0302
pf2 = 4.3835
pf3 = 7.9335