function [xXi,y2]=solve_system(alpha,m,init_coord,approx_period)
%function [xXi]=solve_system(alpha,m)
%
% INPUT:
% alpha           coefs with the system saved
% m               positive integer, number of modes considered in the solution
% init_coord      initial coordinates (must be of compatible size)
% approx_period   approximation of the period of the solution
% OUTPUT:
% x_Xi            Xi vector with m modes solving the Van der Pol equations.

size_vec=alpha.size_vector;

yy=init_coord;
tt1=0;
tt2=1;%approx_period;
t1=linspace(tt1,tt2,2^10*(3*tt2)+1);
ode_func=@(t,y)  general_ode(alpha,t,y)*approx_period*2*pi;
[t2,y2]=ode45(ode_func,t1,yy);

plot(t2,y2);


time=approx_period;
% yy2=y2(t2>=t2(end)-time,:);
% y2=yy2;
% t2=t2(t2>=t2(end)-time);
% t2=t2-t2(1);
x=0*y2;
for i=1:size_vec
    x(:,i)=1/(size(y2,1))*fft(y2(:,i));
end
xBar=cell(size_vec+alpha.size_scalar,1);
xBar{1}=approx_period;%time;
for i=1:size_vec
    xBar{i+alpha.size_scalar}=fftshift([x(1:m+1,i);x(end-m+1:end,i)]);
    xBar{i+alpha.size_scalar}(1:m)=conj(flip(xBar{i+alpha.size_scalar}(m+2:end)));
    xBar{i+alpha.size_scalar}(m+1)=real(xBar{i+alpha.size_scalar}(m+1));
end

size_scal=1;%alpha.size_scalar;
nodes=m;
scal=xBar{1:alpha.size_scalar};
vec=[xBar{1+alpha.size_scalar}];
for i=2:size_vec
    vec=[vec,xBar{i+alpha.size_scalar}];
end
xXi=Xi_vector(scal,vec,size_scal,size_vec,nodes);

return
