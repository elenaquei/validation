function y=F_function(x,alpha,coefs_linear,varargin)
% function y=F_function(x,alpha,coefs_linear,flag)
% 
% INPUT:
% x                  Xi_vector, numerical approximation of the solution;
% alpha              coefs, data of the ODE problem;
% coefs_linear       2-cell with coefficients of scalar part;
% flag (DEFAULT=1)   if flag=1, returns a Xi_vector of the same size of x
%                    else returns a Xi_vector of maximum size.
%
% OUTPUT:
% y                  Xi_vector, value of the function in x.


%%%%%%%%%%THIS FUNCTION NEEDS REWRITING


y=powers(x,alpha,coefs_linear,varargin{:});





%if nargin==3
    % tic;
%    y=powers(x,alpha,coefs_linear);
    % T1=toc;
    % tic;
    % %y2=
    % powers2(x,alpha,coefs_linear); % at the moment, slower than
    % T2=toc;
    % if T1<T2
    %     disp('longer is faster')
    % else
    %     disp(1)
    % end
    % expected.
    % if any(y~=y2)
    %     disp('there is somethign wrong! Do you hear me...')
    % end
%elseif nargin==4
%    y=powers(x,alpha,coefs_linear,flag);
%    y2=powers2(x,alpha,coefs_linear,flag);
%    if any(y~=y2)
%        error('there is somethign wrong2! Do you hear me...')
%    end
%else
%    error('Number of inputs wrong')
%end
return