% example of conitnuation of Hopf

% \dot {x} = - y + x ( mu - x^2 - y^2) 
% \dot {y} =   x + y ( mu - x^2 - y^2)

% rescaling : x -> ax , y -> ay 
% a \dot {x} = - a y + a x mu - a^3 x^3 - a^3 x y^2 
% a \dot {y} =   a x + a y mu - a^3 y x^2 - a^3 y^3


step_size= 10^-4;
number_of_nodes = 4;
number_of_steps = 1;

omega = @(T) 2*pi / T;

f1 = @(a,mu,T,x,y) - a* y + a* x *mu - a.^3* x.^3 - a.^3* x *y.^2 ;
f2 = @(a,mu,T,x,y)   a* x + a* y *mu - a.^3* y* x.^2 - a.^3* y.^3 ;

f1_four = @(a,mu,T,x,y) - a* y + a* x *mu - a^3* conv(x,conv(x,x,'same'),'same') - a^3* conv(x , conv(y,y,'same'),'same') ;
f2_four = @(a,mu,T,x,y)   a* x + a* y *mu - a^3* conv(y,conv(x,x,'same'),'same') - a^3* conv(y , conv(y,y,'same'),'same') ;

dxf1_four = @(a,mu,T,x,y) a*mu *one_four(x) -3*conv(x,x,'same')*a^3-conv(y,y,'same')*a^3;
dxf2_four = @(a,mu,T,x,y) a*one_four(x) - 2*conv(x,y,'same')*a^3;

dyf1_four = @(a,mu,T,x,y) -a *one_four(y) - 2*conv(x,y,'same') *a^3;
dyf2_four =  @(a,mu,T,x,y) mu*a*one_four(x) - conv(x,x,'same') *a^3 - 3*conv(y,y,'same')*a^3;

dmuf1_four =  @(a,mu,T,x,y) x*a;
dmuf2_four =  @(a,mu,T,x,y) y*a;

F1 = @ (a,mu,T,x,y) a * x  - integral_Fourrier_to_F(f1_four(a,mu,T,x,y),omega(T));
F2 = @ (a,mu,T,x,y) a * y - a * one_four(x)- integral_Fourrier_to_F(f2_four(a,mu,T,x,y),omega(T));
F3 = @ (a,mu,T,x,y) integral_Fourrier(f1_four(a,mu,T,x,y),T,omega(T));
F4 = @ (a,mu,T,x,y) integral_Fourrier(f2_four(a,mu,T,x,y),T,omega(T));

F= @(a,mu,T,x,y) [vert(F3(a,mu,T,x,y));vert(F4(a,mu,T,x,y));vert(F1(a,mu,T,x,y));vert(F2(a,mu,T,x,y))];


dxF1 = @ (a,mu,T,x,y) a - integral_Fourrier_to_F(dxf1_four(a,mu,T,x,y),omega(T));
dxF2 = @ (a,mu,T,x,y) - integral_Fourrier_to_F(dxf2_four(a,mu,T,x,y),omega(T));
dxF3 = @ (a,mu,T,x,y) integral_Fourrier_to_F(dxf1_four(a,mu,T,x,y),omega(T));
dxF4 = @ (a,mu,T,x,y) integral_Fourrier_to_F(dxf2_four(a,mu,T,x,y),omega(T));

dyF1 = @ (a,mu,T,x,y) - integral_Fourrier_to_F(dyf1_four(a,mu,T,x,y),omega(T));
dyF2 = @ (a,mu,T,x,y) a - integral_Fourrier_to_F(dyf2_four(a,mu,T,x,y),omega(T));
dyF3 = @ (a,mu,T,x,y) integral_Fourrier_to_F(dyf1_four(a,mu,T,x,y),omega(T));
dyF4 = @ (a,mu,T,x,y) integral_Fourrier_to_F(dyf2_four(a,mu,T,x,y),omega(T));

dmuF1 = @ (a,mu,T,x,y) - integral_Fourrier_to_F(dmuf1_four(a,mu,T,x,y),omega(T));
dmuF2 = @ (a,mu,T,x,y) - integral_Fourrier_to_F(dmuf2_four(a,mu,T,x,y),omega(T));
dmuF3 = @ (a,mu,T,x,y) integral_Fourrier(dmuf1_four(a,mu,T,x,y),T,omega(T));
dmuF4 = @ (a,mu,T,x,y) integral_Fourrier(dmuf2_four(a,mu,T,x,y),T,omega(T));

dTF1 = @ (a,mu,T,x,y) - integral_Fourrier_to_F(f1_four(a,mu,T,x,y),omega(1));
dTF2 = @ (a,mu,T,x,y) - integral_Fourrier_to_F(f2_four(a,mu,T,x,y),omega(1));
dTF3 = @ (a,mu,T,x,y) sum(f1_four(a,mu,T,x,y));
dTF4 = @ (a,mu,T,x,y) sum(f2_four(a,mu,T,x,y));


DF = @ (a,mu,T,x,y) [dmuF3(a,mu,T,x,y)         dTF3(a,mu,T,x,y)           horiz(dxF3(a,mu,T,x,y))        horiz(dyF3(a,mu,T,x,y))
                     dmuF4(a,mu,T,x,y)         dTF4(a,mu,T,x,y)           horiz(dxF4(a,mu,T,x,y))        horiz(dyF4(a,mu,T,x,y))
                     vert(dmuF1(a,mu,T,x,y))   vert(dTF1(a,mu,T,x,y))     conv2mat(dxF1(a,mu,T,x,y))     conv2mat(dyF1(a,mu,T,x,y))
                     vert(dmuF2(a,mu,T,x,y))   vert(dTF2(a,mu,T,x,y))     conv2mat(dxF2(a,mu,T,x,y))     conv2mat(dyF2(a,mu,T,x,y)) ];

% here we should start Newton
% first step: (a, mu, T, x, y) = (0,0,1,sint, cost)
a =0; mu= 0; T=1; y=one_four(zeros(2*number_of_nodes+1,1)); x=0*y; 
x(number_of_nodes)=1/2*1i;x(number_of_nodes+2)=-1/2*1i;

for i=1:number_of_steps
    a=a+step_size;
    [mu,T,x,y] = Newton(@(mu,T,x,y)F(0,mu,T,x,y),@(mu,T,x,y)DF(0,mu,T,x,y),mu,T,x,y);
end

% need
% function integral = integral_Fourrier(c,t,omega)