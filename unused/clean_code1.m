% a nicer code to do the same stuff

% structure of the code:
%  0 ) construct f ( lambda, x) - our vector field
%  1 ) compute Hopf bifurcation and validate it (generalise previous code) -
%  2 ) continuation of the periodic orbit

addpath('./continuation (original)/');

%% construction of the vector field
% we are here using the example presented in the pdf
%
%   dot x1 = -x2 + x1( lambda - x1^2 -x2^2)
%   dot x2 =  x1 + x2( lambda - x1^2 -x2^2)

dim = 2;
deg = [3,3];
non_zero=[4,4];
power = cell(2,1);
power{1} = [0 1
    1 0
    3 0
    1 2];
power{2} = [1 0
    0 1
    2 1
    0 3];
power_lambda = cell(2,1);
power_lambda{1} = [0 1 0 0];
power_lambda{2} = [0 1 0 0];
coef = cell(2,1);
coef{1} = [-1 1 -1 -1];
coef{2} = [1 1 -1 -1];
f = coefficients_polynomial(dim, deg, non_zero, power, power_lambda, coef); % construction of the vector field

%% validation of the Hopf bifurcation
% NOT YET generalised
% for the specific example we know that xH=(0,0) and lambdaH = 0
% analytically ( radius = 0)
x_star =[0,0];
lambda_star = 0;
radius_star =0;
b =1;

%% continuation of the Hopf bifurcation

num_iter = 10; 
step_iter = 10^-4; % step on lambda
n_nodes = 2; % number of Fourier modes

lambda = lambda_star;
x = x_star;
% j structure including the Fourier coefficients of the first
% approximation, a and omega
sin_cos = 1;
p0= [0,1];
[j,f] = strarting_structure(f, x, lambda, sin_cos, b, n_nodes);%(f, sin_cos, b, n_nodes); % updated 

[system, x_xi_vector] = converting_data(j,f);  % converting the data up to now to be compatible with the old code

% call the old code
%NOPE%%%%new_thingy(x_xi_vector.vector,a0, period0,system); % STILL PROBLEMS with multiple scalar unknowns and the scalar equations part!

simplifief_input();




% % j.a = sqrt(step_iter);
% % 
% % for i = 1:num_iter
% %     lambda = lambda + step_iter;
% % %    x = Newton_scalar ( f, lambda, x);   % WORKS - not useful anymore
% %     disp('Launching Newton');
% %     [j, error_out, residual] = Newton_fourier ( f, j);
% %     if error_out ==0
% %         error('Newton did not converge');
% %     elseif error_out ==1
% %         fprintf('Newton converged with residual %f',residual);
% %     elseif error_out ==2
% %         fprintf('Newton converged well with residual  %f',residual);
% %     else
% %         fprintf('Something is weird')
% %     end
% % end
% % 
% % 
% % 

