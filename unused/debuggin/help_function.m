function y = help_function(x)
% x   XI_vector

T= x.scalar(1);
lambda = x.scalar(2);
a = x.scalar(3);
x1 = x.scalar(4);
x2 = x.scalar(5);
y1 = x.vector(1,:);
y2 = x.vector(2,:);
nodes = (length(y1)-1)/2;
Conv = @(x,y) conv(x,y,'same');
der = @(x) 1i*(-nodes:nodes).*x;

first = T *y2 - lambda*T*y1 + 3*T*x1^2*y1+ 3*x1*Conv(y1,y1)*a*T+...
    Conv(y1,Conv(y1,y1))*a^2*T+...
    2*x1*x2*T*y2+a*T*x1*Conv(y2,y2)+...
    T*y1*x2^2+2*x2*a*T*Conv(y2,y1)+a^2*T*Conv(Conv(y1,y2),y2)-der(y1);

second = -der(y2) -T*y1-T*lambda*y2+T*x1^2*y2+2*T*x1*x2*y1+2*x1*T*a*Conv(y1,y2)+x2*T*a*Conv(y1,y1)+...
    T*a^2*Conv(Conv(y1,y1),y2)+3*T*x2^2*y2+3*T*a*x2*Conv(y2,y2)+T*a^2*Conv(Conv(y2,y2),y2);

y = [zeros(1,5),first,second];