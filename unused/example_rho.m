% another example

step_size= 10^-2;
number_of_nodes = 20;
number_of_steps = 1;

inverse_four = @(x) conv2mat(x) \ vert(one_four(x));
hat = @(c) 2*pi*1i*vert(c).*vert(-((length(c)-1)/2):((length(c)-1)/2));
zeroth = @(c) c((length(c)+1)/2);

alpha =@(lambda) lambda;
beta = @(lambda) 1;
g = @(a,b,lambda) - conv(a,conv(a,a,'same'),'same') - conv (a,conv(b,b,'same'),'same');
h = @(a,b,lambda) - conv(a,conv(a,b,'same'),'same') - conv (b,conv(b,b,'same'),'same');

cos=one_four(zeros(2*number_of_nodes+1,1)); sin=0*cos; 
sin(number_of_nodes)=1/2*1i;sin(number_of_nodes+2)=-1/2*1i;

p = @(rho, lambda) conv(conv(rho,rho,'same'),rho,'same');

%g(conv(rho,cos,'same'),conv(rho,sin,'same'),lambda) +...
%    h(conv(rho,cos,'same'),conv(rho,sin,'same'),lambda);
q = @(rho, lambda) 0*rho;
%    conv(h(conv(rho,cos,'same'),conv(rho,sin,'same'),lambda),cos,'same')-...
%    conv(g(conv(rho,cos,'same'),conv(rho,sin,'same'),lambda),sin,'same');


S = @(rho, lambda) lambda * rho - conv(conv(rho,rho,'same'),rho,'same');
%conv(conv(rho,alpha(lambda) * rho + p(rho, lambda),'same'),inverse_four(conv(beta(lambda),rho,'same') - q(rho, lambda)),'same');




F = @(rho,lambda, a) [zeroth(rho) - a;
    vert(hat(rho) - S(rho, lambda))];


DF = @(rho, lambda,a) [0 horiz(one_four(rho));
     vert(-rho),  diag(hat(ones(size(rho))))+conv2mat(-lambda*one_four(rho)+3*conv(rho,rho,'same')) ];


a =0; mu= 0;  rho=one_four(zeros(2*number_of_nodes+1,1)); 

for i=1:number_of_steps
    a=a+step_size;
    [mu,rho] = Newton_rho(@(mu,rho)F(rho,mu,a),@(mu,rho)DF(rho,mu,a),mu,rho);
end