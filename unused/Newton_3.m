function [xBar,yBar,res,DFm,RES]=Newton_3(xBar,alpha,F,DF,maxiter,min_res)
% function [xBar,yBar,res,DFm,RES]=Newton_3(xBar,alpha,F,DF,maxiter,min_res)
% 
% INPUT:
% xBar      Xi_vector, first approximation of the solution in Fourier space;
% alpha     full_problem
% F         inline function
% DF        inline function
% maxiter   positive integer, maximum number of iterations of the Newton
%           method (DEFAULT 20)
% min_res   positive real value, stopping criterion on the residual
%           (DEFAULT 10^-14)
%
% OUTPUT:
% xBar      Xi_vector, final solution;
% yBar      Xi_vector, residual;
% res       positive real value, infinite norm of the residual.
% DFm       last derivative
% RES       sequence of residuals
%
% In the Newton method, the infinite norm on the Xi_vectors is used, that
% is, the biggest value between all is taken, considering scalars and
% elements of the Fourier transforms alike. 
% An error message is triggered if the derivative used is close to
% singularity. 

%global Display % if 1 display the final residual 
global use_intlab
global talkative

if talkative>0
    disp('Entering Newton')
end

if ~exist('maxiter','var') || isempty(maxiter)
    maxiter=10;
end
if ~exist('min_res','var') || isempty(min_res)
    min_res=10^-14;
end

%disp('start')
RES=zeros(maxiter,1);
for i=1:maxiter
    if talkative>1
        fprintf('Iteration %d, time %s\n',i,datestr(now,13));
    end
    yBar=apply(alpha,xBar,0); % OK
    y=Xi_vec2vec(yBar); % working
    x=Xi_vec2vec(xBar);
    res=max(norm(yBar));
    RES(i)=res;
    if talkative>1
        fprintf('Residual %e, time %s\n',res,datestr(now,13));
    end
    %plot(i,res,'*');
    %hold on
    if res<min_res %&& i>1
        break
    end
    
    DFm=derivative_to_matrix(derivative(alpha,xBar,0)); 
    % test
    %max(max(abs(DFm(6:end,6:end) - DF(xBar))))
    %max(abs(y - F(xBar)))
    if ~use_intlab && abs(det((DFm))) > 10e-5
        x=x-DFm\y;
    elseif use_intlab && abs(det(sup(DFm))) > 10e-5
        x=x-DFm\y;
    else
        error('derivative is singular, iteration %d',i)
    end
    xBar=vec2Xi_vec(x,xBar.size_scalar,xBar.size_vector,xBar.nodes);
    xBar = symmetrise(xBar);
end

if i==maxiter && res>min_res
    disp(res)
    return%error('Newton did not converge')
end

if talkative>1
    fprintf('The residual of the Newton method is %d.\n',res)
end
if talkative>0
    fprintf('Exiting Newton\n\n')
end

return