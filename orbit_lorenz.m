function almost_periodic_orbit = orbit_lorenz(hard_mode)
global Display
if nargin == 0
    hard_mode = 0;
end

pho = 28;
sigma = 10;
beta = 8/3;
f1 = @(x) sigma*(x(2) -x(1));
f2 = @(x) x(1)*(pho - x(3)) -x(2);
f3 = @(x) x(1)*x(2) - beta * x(3);
F = @(t,x) [f1(x); f2(x); f3(x)];

y0 = [3.9121    6.6072   14.0230];
[t,y] = ode45(F, [0,300], y0.');

if Display
    plot3(y(:,1),y(:,2),y(:,3))
    hold on
    plot3(y(1,1),y(1,2),y(1,3),'*')
    plot3(y(end,1),y(end,2),y(end,3),'g*')
end

best_index = 5;
best_return = y(best_index, :);
first_point = y(1,:);
i = 5;
if hard_mode
    index_hard = 300;
else
    index_hard = 6;
end
while best_index <index_hard && i < size(y,1)
    local_return = y(i,:);
    if norm(local_return - first_point)<norm(best_return-first_point)
        best_return = local_return;
        best_index = i;
    end
    i = i+1;
end

if Display
    plot3(y(best_index,1),y(best_index,2),y(best_index,3),'k*')
end
figure
plot3(y(1:best_index,1),y(1:best_index,2),y(1:best_index,3))

almost_periodic_orbit = y(1:best_index,:);