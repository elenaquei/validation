function [xBar1,x_dot1,DH1,coefs_linear1] =update_2(xBar0,delta,x_dot0,...
    alpha_coef,alpha_short,coefs_short,coefs_linear,maxiter,min_res)
%function [x1,x_dot1,DH1,coefs_linear1] =update(xBar0,delta,x_dot0,...
%   alpha,coefs_linear,maxiter,min_res)

global use_intlab
temp_intlab=use_intlab;
use_intlab=0;

% compute x_iter
x_dot0_Xi=vec2Xi_vec(x_dot0,xBar0.size_scalar,xBar0.size_vector,xBar0.nodes);
x1=xBar0+ delta* x_dot0_Xi;
%x1.scalar(2)=xBar0.scalar(2)+delta;

% enforce reality of solution
x1=x1.symmetrise;

% add one equation for archlength continuation
coefs_linear0=coefs_linear;
coefs_linear0{1}=[-x_dot0(1:2).';[0,0]];
coefs_linear0{2}(2,:,:)=coefs_linear{2}(1,:,:);
coefs_linear0{2}(1,:,:)=reshape(-x_dot0(3:end),2*xBar0.nodes+1,xBar0.size_vector).';
coefs_linear0{3}(1) = abs(Xi_vec2vec(x1).'*x_dot0);
coefs_linear0{3}(2)=0;
coefs_linear0{3}(2)=real(-coefs_linear0{1}(2,:)*x1.scalar.'-...
    sum(sum(squeeze(coefs_linear0{2}(2,:,:)).*x1.vector,1)));


% Newton
%x1short=x1;
%x1short.size_scalar=1;
%x1short.scalar=x1short.scalar(1);
%alpha_short.value{1}(5)=-(x1.scalar(2)+delta)/(2*pi);
%[xBar1]=Newton_Xi(x1short,alpha_short,coefs_short,maxiter,min_res);
[xBar1]=Newton_Xi(x1,alpha_coef,coefs_linear0,maxiter,min_res);
% impose symmetrisation (?)
xBar1=xBar1.symmetrise;
%xBar1.size_scalar=2;
%xBar1.scalar(2)=x1.scalar(2)+delta;


DH1=Function_derivative(xBar1,alpha_coef,coefs_linear,0);
% flag 0 to have non-square output

x_dot1=null(DH1);


%x_dot1=x_dot1 *conj(el);

% symmetric = @(x) (x+conj(x(end:-1:1)))/2;
% 
% x_dot1(1:alpha_coef.size_scalar)= abs(x_dot1(1:alpha_coef.size_scalar));
% for ii=1:alpha_coef.size_vector
%     indeces= alpha_coef.size_scalar+ (ii-1)*(2*xBar1.nodes+1) +...
%     1:(2*xBar1.nodes+1);
%     el=x_dot1( alpha_coef.size_scalar+ (ii-1)*(2*xBar1.nodes+1) + xBar1.nodes+1);
%     x_dot1(indeces) = (conj(el)/abs(el)*(x_dot1(indeces));
% end
% x_dot1=x_dot1/norm(x_dot1);

if dot(x_dot1,x_dot0)<0
    x_dot1=-x_dot1;
end
theta_best=0;
for theta=0:0.0001:2*pi 
    if norm(x_dot0-exp(1i*theta)*x_dot1)<norm(x_dot0-exp(1i*theta_best)*x_dot1)...
            && dot(exp(1i*theta)*x_dot1,x_dot0)>0
        theta_best=theta;
    end
end
x_dot1=exp(1i*theta_best)*x_dot1;

% el=x_dot1( alpha_coef.size_scalar+ xBar0.nodes+1);
% x_dot1=x_dot1 *conj(el)/abs(el);
% if dot(x_dot1,x_dot0)<0
%     x_dot1=-x_dot1;
% end

if max(abs(DH1*x_dot1)) > 10^-10
    error('inversion not cool')
end


if norm(x_dot0-x_dot1) > 100*delta
    warning('Norm derivatives in x0 and x1 quite big')
elseif any(norm(xBar0-xBar1) < delta/100)
    warning('Distance between consecutive points quite small, %e',norm(xBar0-xBar1))
end



coefs_linear1=coefs_linear;
coefs_linear1{1}=[-x_dot1(1:2).';[0,0]];
coefs_linear1{2}(2,:,:)=coefs_linear1{2}(1,:,:);
coefs_linear1{2}(1,:,:)=reshape(-x_dot1(3:end),2*xBar0.nodes+1,xBar0.size_vector).';
coefs_linear1{3}(1) = Xi_vec2vec(xBar1).'*x_dot1;
coefs_linear1{3}(2)=real(-coefs_linear1{1}(2,:)*xBar1.scalar.'-...
    sum(sum(squeeze(coefs_linear1{2}(2,:,:)).*xBar1.vector,1)));



DH1=[-x_dot1.';DH1];


use_intlab=temp_intlab;
return
