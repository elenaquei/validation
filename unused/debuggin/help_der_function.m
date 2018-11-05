function DF = help_der_function(x)
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
der = 1i*(-nodes:nodes);
one_four = 0*der;
one_four(nodes+1)=1;

D1F1 = -T*lambda*one_four + 3*T*x1^2*one_four+6*T*a*x1*y1+3*T*a^2*Conv(y1,y1)+...
    T*x2^2*one_four+2*T*x2*a*y2+T*a^2*Conv(y2,y2);

D2F1 = (-T+2*T*x1*x2)*one_four +2*T*x1*a*y2 +2*T*a*x2*y1+2*T*a^2*Conv(y1,y2);

D1F2 = T*one_four +2*T*x2*a*y1+2*T*a^2*Conv(y1,y2)+2*x1*x2*T*one_four+T*a*x1*y2;

D2F2 = (- lambda*T+T*x1^2)*one_four+2*T*x1*a*y1+T*a^2*Conv(y1,y1)+3*T*x2^2+6*T*a*x2*y2+...
    3*T*a^2*Conv(y2,y2);

pad =zeros(1,nodes);

DF = [-1i*diag(-nodes:nodes)+toeplitz([D1F1(nodes+1:-1:1),pad],[D1F1(nodes+1:end),pad])   toeplitz([D2F1(nodes+1:-1:1),pad],[D2F1(nodes+1:end),pad])
toeplitz([D1F2(nodes+1:-1:1),pad],[D1F2(nodes+1:end),pad])   -1i*diag(-nodes:nodes)+toeplitz([D2F2(nodes+1:-1:1),pad],[D2F2(nodes+1:end),pad])];

