% TEMP_MAIN HOPF
%
% a temporary main with minimum functionality to start somewhere to built
% the code

% STRUCTURE:
% -) set alpha, the full problem
% -) find a numerical periodic solution
% -) apply Newton_Xi to get a better one
% -) compute the radii polynomials
% -) continue

%% full problem
% here we define frist and solve the Hopf problem

%   dot x1 = -x2 + x1( lambda - x1^2 -x2^2)
%   dot x2 =  x1 + x2( lambda - x1^2 -x2^2)

global azabaza
global nu
global use_intlab 
global talkative 
global RAD_MAX
talkative = 1;
use_intlab = 0;
nu = 1.1;
RAD_MAX = 10^-4;

if isempty(azabaza)
    addpath(genpath('../'))
    startintlab
    azabaza =1;
end

n_nodes= 3;

dim_vec = 2;
dim_scal = 1;
n_eq = dim_vec;
non_zero=[4,4];
power = cell(2,1);
power{1}{1} = [0 1]';
power{1}{2} = [1 0]';
power{1}{3} = [3 0]';
power{1}{4} = [1 2]';

power{2}{1} = [1 0]';
power{2}{2} = [0 1]';
power{2}{3} = [2 1]';
power{2}{4} = [0 3]';

power_lambda = cell(2,1);
power_lambda{1} = [0 1 0 0];
power_lambda{2} = [0 1 0 0];
coef = cell(2,1);
coef{1} = [-1 1 -1 -1];
coef{2} = [1 1 -1 -1];


% %% DUM TEST 
% % f(x) = (x1 +LAMBDA * x1 * x2  ;x2 + x2 * x2)
% 
% dim_vec = 2;
% dim_scal = 1;
% n_eq = dim_vec;
% non_zero=[2,2];
% power = cell(2,1);
% power{1}{2} = [1 1]';
% power{1}{1} = [1,0]';
% 
% power{2}{1} = [0 1]';
% power{2}{2} = [0 2]';
% 
% power_lambda = cell(2,1);
% power_lambda{1} = [0, 1];
% power_lambda{2} = [0, 0];
% coef = cell(2,1);
% coef{1} = [1,1];
% coef{2} = [1,1];
% 
f = polynomial_coefs(dim_scal, dim_vec, n_eq, non_zero, coef, power_lambda, power);
F = Taylor_series_Hopf(f,n_nodes); % debugged

%% solution

x_star =[0,0];
lambda_star = 0;
T_star = 1;%/(2*pi);  
a_star = 0;
sin_four = zeros(1,2*n_nodes+1);
cos_four = sin_four;
cos_four(n_nodes) = 1/2; cos_four(n_nodes+2) = 1/2;
sin_four(n_nodes) = 1i/2; sin_four(n_nodes+2) = -1i/2;

y = [sin_four; cos_four];

added_err = 10^-3;

x =  x_star;%added_err*rand
lambda = lambda_star; %added_err*rand
T = T_star; %added_err*rand
a = added_err+ a_star;

sol_star = Xi_vector([T_star, lambda_star, a_star, x_star],y);
sol = Xi_vector([T, lambda, a, x],y);

% residual = apply(F,sol_star);
% DF =  derivative(F,sol_star);    % scalar derivative perfectly right
% vector derivative matches in dimensions, but something probably off
% DF_mat = derivative_to_matrix(DF);

p1 = Xi_vec2vec(sol_star);
lin_coef = cell(3,1);

lin_coef{1} = 0*sol_star.scalar;
lin_coef{2} = shiftdim(sol_star.vector,-1);
lin_coef{2}(1,1,:) = sum(sin_four) + 0*sin_four;
lin_coef{2}(1,2,:) = sum(cos_four) + 0*sin_four;
lin_coef{3} = -1;
F.scalar_equations = change_lin_coef(F.scalar_equations,lin_coef,2);

lin_coef{1} = 0* sol_star.scalar;
lin_coef{2} = 1+0*shiftdim(sol_star.vector,-1);
der1 = (1i*(-n_nodes:n_nodes).*y(1,:));
der2 = (1i*(-n_nodes:n_nodes).*y(2,:));
lin_coef{2}(1,1,:) = sum(der1);
lin_coef{2}(1,2,:) = sum(der2);
lin_coef{3} = -sum(p1(6:end));
F.scalar_equations = change_lin_coef(F.scalar_equations,lin_coef,3);

F.scalar_equations.linear_coef{3}(end) = -a;

%f_inline = @(x) help_function(x);
%Df_inline = @(x) help_der_function(x);
F.vector_field.value{1}(end) =-1;
F.vector_field.value{2}(end) =-1;

[sol2,yBar,res,DFm,RES] =Newton_2(sol,F,30,10^-6); % IT WOOORKSSSSS !!!! :D :D :D 
sol2.scalar(3) = sol2.scalar(3);% + added_err^4*rand;
% sol2 = symmetrise(sol2);
%res = norm(apply(F,sol2));
% %% DUM TEST 2 - for DF - seems to be working, or at least this test
% doesn't detect any error.
% % the scalar part is working with no trouble, so I need a test about just
% % the vector part 
% % and keeping things simple:
% % f(x1,x2) = (lambda * x1+x2+x2^2
% %              x1 * x2)
% 
% dim_vec = 2;
% dim_scal = 1;
% n_eq = dim_vec;
% non_zero=[3,1];
% power = cell(2,1);
% power{1}{1} = [1,0]';
% power{1}{2} = [0 1]';
% power{1}{3} = [0 2]';
% 
% power{2}{1} = [1 1]';
% 
% power_lambda = cell(2,1);
% power_lambda{1} = [1,0,0];
% power_lambda{2} = 0;
% coef = cell(2,1);
% coef{1} = [1,1,1];
% coef{2} = [1];
% dot = power;
% dot{1}{1} = [0,0]';
% dot{1}{2} = [0 0]';
% dot{1}{3} = [0 0]';
% 
% dot{2}{1} = [0 0]';
% 
% sol = Xi_vector(1,y);
% f = polynomial_coefs(dim_scal, dim_vec, n_eq, non_zero, coef, power_lambda, power,dot);
% [Dlambda, Dx_diag, Dx_vec] = compute_derivative(f,sol); % right
% 
% poly = polynomial_coefs(1, 2, 0, ...
%                 [], cell(0,0),cell(0,0),cell(0,0), cell(0,0), cell(0,0));
% lin_coef = cell(3,1);
% 
% lin_coef{1} = sol.scalar;
% lin_coef{2} = shiftdim(sol.vector,-1);
% lin_coef{3} = -1;
% 
% g = scalar_eq(1, 0, 1, 2, lin_coef, poly);
% F = full_problem(g,f);
% D=derivative(F,sol);

%% validation
DF =  derivative(F,sol2,0);    
DF_mat = derivative_to_matrix(DF);
A  = inv(DF_mat);

Y_vector = Y_bound_new(A,sol2,F);
Z0_vector=Z0_bound(DF_mat,A,sol2);
Z1_vector=Z1_bound_new(A,sol2,F);
Z2_vector= Z2_bound_new(A,sol2,F);
[Imin,Imax]=find_negative(Z2_vector,Z1_vector,Z0_vector,Y_vector);


% DF_non_square = DF_mat([1,2,4:size(DF_mat,1)],:);
% v = kernel(DF_non_square);
% DF_new = DF_mat;
% DF_new(3,:) = v;

% %% continuation step
% A2 = A;
% DF_mat2 = DF_mat;
% F2 = F;
% a = 10*added_err +a;
% F3 = F2;
% F3.scalar_equations.linear_coef{3}(end) = -a;
% sol3 = sol2;
% sol3.scalar(3) = a;
% [sol3,yBar,res,DFm,RES] =Newton_2(sol3,F3,30,10^-6);
% 
% DF3 =  derivative(F3,sol3,0);    
% DF_mat3 = derivative_to_matrix(DF3);
% A3  = inv(DF_mat3);
% 
% Y_vector3 = Y_bound_new(A,sol3,F3);
% Z0_vector3=Z0_bound(DF_mat,A,sol3);
% Z1_vector3=Z1_bound_new(A,sol3,F3);
% Z2_vector3= Z2_bound_new(A,sol3,F3);
% [Imin3,Imax3]=find_negative(Z2_vector3,Z1_vector3,Z0_vector3,Y_vector3);
% 
% %% continuated validation
% Ycont = Y_bound_cont_new(A2,A3,sol2,sol3,F2,F3,Y_vector);
% Z0cont = Z0_bound_cont(DF_mat2,DF_mat3,A2,A3,sol2);
% Z1cont = Z1_bound_cont_new(A2,A3,sol2,sol3,F2,F3,DF_mat2 - DF_mat3,Z1_vector);
% Z2cont = Z2_bound_cont_new(A2,A3,sol2,sol3,F2);

% IT WORKS

%% serious continuation
% starting from sol2, going until -sol2

max_iter = 60;
h = 0.5*10^-2;
a0 = -0.1;
a_old = a0;
a_new = a_old;
Fold = F;
Fold.scalar_equations.linear_coef{3}(end) = -a0;

solold =Newton_2(sol,Fold,30,10^-10);
sol0=solold;
DFold =  derivative(Fold,solold,0);    
DF_matold = derivative_to_matrix(DFold);
Aold  = inv(DF_matold);


Y_old = Y_bound_new(Aold,solold,Fold);
Z0_old=Z0_bound(DF_matold,Aold,solold);
Z1_old=Z1_bound_new(Aold,solold,Fold);
Z2_old= Z2_bound_new(Aold,solold,Fold);
[Imin3,Imax3]=find_negative(Z2_old,Z1_old,Z0_old,Y_old);


size_x = solold.size_scalar+ solold.size_vector;
previous_iter.Y=[];
previous_iter.Z1=[];
delta_vec=zeros(1,max_iter);
Imin_vec=zeros(1,max_iter);
Imax_vec=zeros(1,max_iter);
xDelta= zeros(1,max_iter);
Y0_vec=zeros(size_x,max_iter);
Z0_vec=zeros(size_x,max_iter);
Z1_vec=zeros(size_x,max_iter);
Z2_vec=zeros(size_x,max_iter);

for iter = 1:max_iter
    if abs(a_new - a0) > 2*abs(a0)
    %   break
    end
    a_new = a_old + h;
    if talkative>0
        fprintf('\nThe amplitude is     %f\n',a_new)
    end
    if talkative>2
        fprintf('Iteration number %d started\n', iter)
    end
    Fnew = Fold;
    Fnew.scalar_equations.linear_coef{3}(end) = -a_new;
    solnew =Newton_2(solold,Fnew,30,10^-10);
    DFnew =  derivative(Fnew,solnew,0);
    DF_matnew = derivative_to_matrix(DFnew);
    Anew  = inv(DF_matnew);
    
    % normal bounds
    if talkative>1
        disp('Started the computation of the bounds')
    end
    [flag,Imin,Imax,previous_iter,Yvector,Z0vector,Z1vector,Z2vector,new_step] = radii_polynomials_cont_new(solold,solnew,DF_matold,DF_matnew,...
    Fold,Fnew,previous_iter,Aold,Anew);
    %h = h*new_step;
    if flag==0
        break
    end
%     Y_new = Y_bound_new(Anew,solnew,Fnew);
%     Z0_new=Z0_bound(DF_matnew,Anew,solnew);
%     Z1_new=Z1_bound_new(Anew,solnew,Fnew);
%     Z2_new= Z2_bound_new(Anew,solnew,Fnew);
%     find_negative(Z2_new,Z1_new,Z0_new,Y_new);
% %     
% %     % continuation bounds
%     Ycont = Y_bound_cont_new(Aold,Anew,solold,solnew,Fold,Fnew);
%     Z0cont = Z0_bound_cont(DF_matold,DF_matnew,Aold,Anew,solold);
%     Z1cont = Z1_bound_cont_new(Aold,Anew,solold,solnew,Fold,Fnew,DF_matold - DF_matnew);
%     Z2cont = Z2_bound_cont_new(Aold,Anew,solold,solnew,Fold);
%     [Imin,Imax]=find_negative(Z2cont,Z1cont,Z0cont,Ycont);
    if talkative>0
        disp('The bound was validated successfully.')
    end
    if talkative>1
        fprintf('with bounds [%e , %e]',Imin, Imax);
    end
    
    % storage
    delta_vec(iter)=solnew.scalar(2);
    Imin_vec(iter)=Imin;
    Imax_vec(iter)=Imax;
    xDelta(iter) = norm_Xi_vector(solold-solnew,nu);
    Y0_vec(:,iter)=Yvector;
    Z0_vec(:,iter)=Z0vector;
    Z1_vec(:,iter)=Z1vector;
    Z2_vec(:,iter)=Z2vector;
    
    % update
    a_old = a_new;
    Fold = Fnew;
    solold = solnew;
    Aold = Anew;
    DF_matold = DF_matnew;
    
    
end

alpha = Fnew;
%alpha.scalar_e
alpha.scalar_equations = remove_lin_coef(alpha.scalar_equations,3);
DH0 = derivative_to_matrix(derivative(alpha,solold,0));
x_dot0=kernel(DH0);
[solnew,x_dot1,DF_matnew,alphanew] =update_new(solold,h,x_dot0,...
    alpha,10,10^-7); % WORKS
