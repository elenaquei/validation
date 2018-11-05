% MAIN_ROSSLER
% as of 31 August 2017: it runs
% as of 4th July 2018: it doesn't - remodernizing the code
%
% continuing the Rossler periodic orbit outside of the origin
%
% STRUCTURE:
% -) set alpha, the full problem
% -) find a numerical periodic solution
% -) apply Newton_Xi to get a better one
% -) compute the radii polynomials
% -) continue

%% full problem
% here we define frist and solve the Hopf problem

% dot x = A*(x-y)-y*z+w;
% dot y = -B*y+x*z
% dot z = -C*z+D*x+x*y
% dot w = -E*(x+y)];
% 
% A parameter 
% B = 4;
% C = B;
% D = 0.04;
% E = 1.4;
% 
% proven Hopf bifurcation at A = 0, (x,y,z,w) = (0,0,0,0)

global nu
global use_intlab 
%global talkative 
global RAD_MAX
%talkative = 3;
use_intlab = 0;
nu = 1.1;
RAD_MAX = 10^-4;

try 
    intval(1);
catch
    addpath(genpath('../'));
    startintlab;
end
% if isempty(azabaza)
%     addpath(genpath('../'))
%     startintlab
%     azabaza =1;
% end

hmin = 10^-4;

B = 4;
C = B;
D = 0.04;
E = 1.4;

dim_vec = 4;
dim_scal = 1;
n_eq = dim_vec;

% dot x = A*x  -A*y -y*z +w
% dot y = -B*y +x*z
% dot z = -C*z +D*x +x*y
% dot w = -E*x +Ey 

string_rossler = '- dot x1 + A x1 - A x2 - x2 x3 + x4 \n - dot x2 -B x2 + x2 x3 \n - dot x3 - C x3 + D x1 + x1 x2 \n - dot x4 - E x1 + E x2'; % Rossler ( weird one)
string_rossler_vars = strrep(string_rossler, 'B' , num2str(B)); % plugging in the numerical value for B
string_rossler_vars = strrep(string_rossler_vars, 'C' , num2str(C)); % plugging in C
string_rossler_vars = strrep(string_rossler_vars, 'D' , num2str(D)); % plugging in D
string_rossler_vars = strrep(string_rossler_vars, 'E' , num2str(E)); % plugging in E

string_rossler_cont = strrep(string_rossler_vars, 'A', 'l1'); % setting A as the parameter for the Hopf bifurcation




% some elements useful for the computation and the validation
n_nodes = 5; % number of Fourier nodes used: small, since near the Hopf bifurcation is a circle
n_iter = 5; % number of iterations
step_size = 10^-4; % initial step size (then adapted along the validation
s = 'hopf_Rossler'; % where the solutions are stored

%vectorfield = strrep(string_lorenz84_cont, 'l1' , '');%'-dot x1 - x2 + l1 x1 - x1 ^ 3 - x1  x2 ^ 2\n- dot x2 + x1 + l1 x2 - x1 ^ 2 x2 - x2 ^ 3';
%vectorfield = strrep(string_lorenz84_cont, 'l2' , 'l1');
% string defining the vector field of the Hopf normal form 

f_ros = from_string_to_polynomial_coef(string_rossler_cont); % transformation into a vectorfield that can be used

% definition of the solution
dummy = load('hopf_in_Rossler');

% validated values converted to doubles 
eigenval = conj(mid(dummy.eigenval));
eigenvec = conj(mid(dummy.eigenvec));
lambda_star = dummy.lambda_star;
x_star = dummy.x_star;
sign_FLC = sign(mid(dummy.l1));
%eigenvec = eigenvec / norm(eigenvec);

% starting the continuation
[s, last_sol] = continuation_Hopf( lambda_star, x_star, f_ros, n_nodes, n_iter, step_size, s, eigenvec, eigenval, sign_FLC);





return






% non_zero=[4,2,3,2];
% power = cell(4,1);
% power{1}{1} = [1 0 0 0]';
% power{1}{2} = [0 1 0 0]';
% power{1}{3} = [0 1 1 0]';
% power{1}{4} = [0 0 0 1]';
% 
% power{2}{1} = [0 1 0 0]';
% power{2}{2} = [1 0 1 0]';
% 
% power{3}{1} = [0 0 1 0]';
% power{3}{2} = [1 0 0 0]';
% power{3}{3} = [1 1 0 0]';
% 
% power{4}{1} = [1 0 0 0]';
% power{4}{2} = [0 1 0 0]';
% 
% 
% power_lambda = cell(4,1);
% power_lambda{1} = [1 1 0 0];
% power_lambda{2} = [0 0 0 0];
% power_lambda{3} = [0 0 0 0];
% power_lambda{4} = [0 0 0 0];
% 
% coef = cell(4,1);
% coef{1} = [1 1 -1 1];
% coef{2} = [-B 1];
% coef{3} = [-C D 1];
% coef{4} = [-E E];
% 
% f = polynomial_coefs(dim_scal, dim_vec, n_eq, non_zero, coef, power_lambda, power);
% 
f = from_string_to_polynomial_coef(vectorfield);
f = set_zero_derivative(f);
F = Taylor_series_Hopf(f,n_nodes); 

%% solution

x_star =[0,0,0,0];
lambda_star = 0;
T_star = 1;%/(2*pi);  
a_star = 0;
sin_four = zeros(1,2*n_nodes+1);
cos_four = sin_four;
cos_four(n_nodes) = 1/2; cos_four(n_nodes+2) = 1/2;
sin_four(n_nodes) = 1i/2; sin_four(n_nodes+2) = -1i/2;

y = [sin_four; cos_four; sin_four; cos_four]; %%%% LOOK FOR EIGS AND EIGENVEC
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
F.scalar_equations = change_lin_coef(F.scalar_equations,lin_coef,1);

lin_coef{1} = 0* sol_star.scalar;
lin_coef{2} = 1+0*shiftdim(sol_star.vector,-1);
der1 = (1i*(-n_nodes:n_nodes).*y(1,:));
der2 = (1i*(-n_nodes:n_nodes).*y(2,:));
lin_coef{2}(1,1,:) = sum(der1);
lin_coef{2}(1,2,:) = sum(der2);
lin_coef{3} = -sum(p1(6:end));
F.scalar_equations = change_lin_coef(F.scalar_equations,lin_coef,2);

F.scalar_equations.linear_coef{3}(end) = -a;

%f_inline = @(x) help_function(x);
%Df_inline = @(x) help_der_function(x);
%F.vector_field.value{1}(end) =-1;
%F.vector_field.value{2}(end) =-1;

[sol2,yBar,res,DFm,RES] =Newton_2(sol,F,30,10^-6); % IT WOOORKSSSSS !!!! :D :D :D 
sol2.scalar(3) = sol2.scalar(3);% + added_err^4*rand;

%% validation
DF =  derivative(F,sol2,0);    
DF_mat = derivative_to_matrix(DF);
A  = inv(DF_mat);

Y_vector = Y_bound_new(A,sol2,F);
Z0_vector=Z0_bound(DF_mat,A,sol2);
Z1_vector=Z1_bound_new(A,sol2,F);
Z2_vector= Z2_bound_new(A,sol2,F);
[Imin,Imax]=find_negative(Z2_vector,Z1_vector,Z0_vector,Y_vector);


%% serious continuation
% starting from sol2, going until -sol2

max_iter = 50;
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
xDelta= zeros(size_x,max_iter);
x_iter= zeros(size_x,max_iter);
Y0_vec=zeros(size_x,max_iter);
Z0_vec=zeros(size_x,max_iter);
Z1_vec=zeros(size_x,max_iter);
Z2_vec=zeros(size_x,max_iter);


alpha = Fold;
%alpha.scalar_e
alpha.scalar_equations = remove_lin_coef(alpha.scalar_equations,3);

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
    
%     DH0 = derivative_to_matrix(derivative(alpha,solold,0));
%     x_dot0=kernel(DH0);
%     [solnew,x_dot1,DF_matnew,Fnew] =update_new(solold,h,x_dot0,...
%         alpha,10,10^-7);
    
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
    
     Y_new = Y_bound_new(Anew,solnew,Fnew);
     Z0_new=Z0_bound(DF_matnew,Anew,solnew);
     Z1_new=Z1_bound_new(Anew,solnew,Fnew);
     Z2_new= Z2_bound_new(Anew,solnew,Fnew);
     find_negative(Z2_new,Z1_new,Z0_new,Y_new);
    
    [flag,Imin,Imax,previous_iter,Yvector,Z0vector,Z1vector,Z2vector,new_step] = radii_polynomials_cont_new(solold,solnew,DF_matold,DF_matnew,...
    Fold,Fnew,previous_iter,Aold,Anew);
    %h = h*new_step;
    if flag==0
        h = h/2;
        if h <hmin
            error('Step size decreased too much')
        else
            continue
        end
    end
    
    if talkative>0
        disp('The bound was validated successfully.')
    end
    if talkative>1
        fprintf('with bounds [%e , %e]',Imin, Imax);
    end
    if talkative>2
        fprintf('Iteration %n out of %n',iter, max_iter);
    end
    % storage
    delta_vec(iter)=solnew.scalar(2);
    Imin_vec(iter)=Imin;
    Imax_vec(iter)=Imax;
    xDelta(:,iter) = norm(solold-solnew);
    x_iter(:,iter) = norm(solnew);
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



% cool pictures in matlab
% axes('FontSize',20)
% plot(x_iter(2,1:20),x_iter(3,1:20),'LineWidth',2)
% hold on
% plot(x_iter(2,20:iter),-x_iter(3,20:iter),'LineWidth',2)
% xlabel('p','Interpreter','Latex','FontSize',20)
% ylabel('a','Interpreter','Latex','FontSize',20)
% title('Parameter and amplitude','Interpreter','Latex','FontSize',20)



 % WORKS