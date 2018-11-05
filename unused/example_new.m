% some other ways of approaching the problem

lambda_star =0;
step_size = 0.001;
number_of_nodes = 3;


Conv = @(x,y) conv(x,y,'same');

F = @(x,y,a,l,o) [ (( l - 3*x(1)^2 -x(2)^2)*y(1)+ (-1-2*x(1)*x(2))*y(2) +...
    (-6* x(1))*y(1)^2 * a+ 2* (-2*x(2))*y(1)*y(2) * a + (-2*x(1))*y(2)^2*a + ...
    (-6)*y(1)^3*a^2+ 3*(-2) *y(1)*y(2)^2*a^2);
    ((1-2*x(1)*x(2))*y(1) + (l - 3*x(2)^2 - x(1)^2)*y(2) + (-2*x(2))*y(1)^2*a+...
    2*(-2*x(1))*y(1)*y(2)*a +(-6*x(2))*y(2)^2*a + 3*(-2)*y(1)^2*y(2)*a^2+(-6)*y(2)^3*a^2)];

% assuming x(1,:) is a Fourier series and so on for x(2,:), y(1,:)....
F_four_1 = @(x,y,a,l,o)  vert( (l- 3*x(1)^2 -x(2)^2)*y(1,:)+ (-1 -2*x(1)*x(2) )*y(2,:) +...
    (-6)*x(1)*Conv(y(1,:),y(1,:)) * a+ 2*(-2*x(2)) *Conv(y(1,:),y(2,:)) * a + (-2*x(1))*Conv(y(2,:),y(2,:))*a + ...
    (-6)*Conv(Conv(y(1,:),y(1,:)),y(1,:))*a^2 + 3*(-2) *Conv(Conv(y(1,:),y(2,:)),y(2,:))*a^2 );

F_four_2 = @(x,y,a,l,o) vert( ( 1-2*x(1)*x(2))*y(1,:) + (l*x(1) - 3*x(2)*x(2) - x(1)*x(1))*y(2,:) + (-2*x(2))*Conv(y(1,:),y(1,:))*a+...
    2*(-2*x(1))*Conv(y(1,:),y(2,:))*a +(-6*x(2))*Conv(y(2,:),y(2,:))*a + 3*(-2)*Conv(Conv(y(1,:),y(1,:)),y(2,:))*a^2+(-6)*Conv(Conv(y(2,:),y(2,:)),y(2,:))*a^2 );

% F_four_1 = @(x,y,a,l,o) [ vert(( Conv( l*one_four(x(1,:))- 3*Conv(x(1,:),x(1,:)) -Conv(x(2,:),x(2,:)),y(1,:))+ Conv(-one_four(x(1,:))-2*Conv(x(1,:),x(2,:)),y(2,:)) +...
%     (-6)*Conv(x(1,:),Conv(y(1,:),y(1,:))) * a+ 2* Conv(Conv(-2*x(2,:),y(1,:)),y(2,:)) * a + Conv(Conv(-2*x(1,:),y(2,:)),y(2,:))*a + ...
%     (-6)*Conv(Conv(y(1,:),y(1,:)),y(1,:))*a^2+3*(-2) *Conv(Conv(y(1,:),y(2,:)),y(2,:))*a^2))];
% 
% F_four_2 = @(x,y,a,l,o) vert((Conv(1*one_four(x(1,:))-2*Conv(x(1,:),x(2,:)),y(1,:)) + Conv(l*one_four(x(1,:)) - 3*Conv(x(2,:),x(2,:)) - Conv(x(1,:),x(1,:)),y(2,:)) + (-2)*Conv(Conv(x(2),y(1)),y(1,:))*a+...
%     2*Conv(Conv(-2*x(1,:),y(1,:)),y(2,:))*a +Conv(Conv(-6*x(2,:),y(2,:)),y(2,:))*a + 3*(-2)*Conv(Conv(y(1,:),y(1,:)),y(2,:))*a^2+(-6)*Conv(Conv(y(2,:),y(2,:)),y(2,:))*a^2));
F_four = @(x,y,a,l,o) [F_four_1(x,y,a,l,o);
    F_four_2(x,y,a,l,o)];


%[ vert(( Conv( l+ 3*Conv(x(1,:),x(1,:)) -Conv(x(2,:),x(2,:)),y(1,:))+ Conv(-1-2*Conv(x(1,:),x(2,:)),y(2,:)) +...
%    (-6)*Conv(y(1,:),y(1,:)) * a+ 2* Conv(Conv(-2*x(2,:),y(1,:)),y(2,:)) * a + Conv(Conv(-2*x(1),y(2)),y(2,:))*a + ...
%    (-6)*Conv(Conv(y(1,:),y(1,:)),y(1,:))*a^2+3*(-2) *Conv(Conv(y(1,:),y(2,:)),y(2,:))*a^2));
%    vert((Conv(1-2*Conv(x(1,:),x(2,:)),y(1,:)) + Conv(l - 3*Conv(x(2,:),x(2,:)) - Conv(x(1,:),x(1,:)),y(2,:)) + (-2)*Conv(Conv(x(2),y(1)),y(1,:))*a+...
%    2*Conv(Conv(-2*x(1,:),y(1,:)),y(2,:))*a +Conv(Conv(-6*x(2,:),y(2,:)),y(2,:))*a + 3*(-2)*Conv(Conv(y(1,:),y(1,:)),y(2,:))*a^2+(-6)*Conv(Conv(y(2,:),y(2,:)),y(2,:))*a^2))];



f= @(x,l) [ ( -x(2) + x(1)*(l-x(1)^2 -x(2)^2));
    x(1) + x(2) *(l-x(1)^2-x(2)^2)];

Df = @(x,l) [ (l - 3*x(1)^2 - x(2)^2)    (-1-2*x(1)*x(2))
    (1-2*x(1)*x(2))      (l-3*x(2)^2 - x(1)^2)];


lambda = lambda_star + step_size;

x_star = Newton_small ( @(x) f(x,lambda), @(x) Df(x,lambda), [0;0]);

all_eigs=eig(Df(x_star,lambda_star));
[~,index_b]=min(real(all_eigs));
b=abs(imag(all_eigs(index_b)));
omega = 1 /b;
% converging

% now the problem gets to be
% y_dot = F_four          ODE
% y(1) = x_star(1)        phase condition
% (F_four)_0 = 0          periodicity condition

line_zeroth = @(vec) vec( (size(vec,1)+1)/2,:,:,:);

hat = @(c,omega) omega*1i*vert(c).*vert(-((length(c)-1)/2):((length(c)-1)/2));

big_problem= @(y,a,o) [ vert([vert(hat(y(1,:),o));vert(hat(y(2,:),o))])  - vert(F_four(x_star,y,a,lambda,o))
    sum(y(1,:)) - x_star(1);
    line_zeroth(vert(F_four_1(x_star,y,a,lambda,o)))];

y_cos_F=0*one_four(zeros(2*number_of_nodes+1,1)); 
y_cos_F(number_of_nodes) = 1/2; y_cos_F(number_of_nodes+2) = 1/2;
y_sin_F=0*y_cos_F; 
y_sin_F(number_of_nodes)=-1/2*1i;y_sin_F(number_of_nodes+2)=1/2*1i;

y = [horiz(y_sin_F);horiz(y_cos_F)];

big_norm=norm(big_problem(y,0,omega));
sol = big_problem(y,0,omega);

stupid_problem = @(y,o)  [vert(hat(y(1,:),o)) + b* vert(y(2,:)) ; vert(hat(y(2,:),o)) - b* vert(y(1,:)) ];

stupid_norm = norm ( stupid_problem(y,omega));
