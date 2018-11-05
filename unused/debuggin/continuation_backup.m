function [s, x_n] = continuation_backup ( x0, F_not_square, n_iter, h, x_dot_0,s, min_res_N, bool_Hopf)
% function [s, x_n, a_lot_of_stuff] = continuation ( x0, F_not_square, n_iter, h, x_dot_0,s,min_res_N, bool_Hopf)
% WOEKS BUT JUST FIRST ROUND
%
% INPUT
% x0             Xi_vector, solution at one point
% F_not_square   full_problem, underdetermined
% n_iter         integer, number of interations (DEFAULT: 10)
% h              real, initial step size (DFAULT: 10^-5)
% x_dot_0        complex vector of length M+N*(nodes*2+1), kernel of
%                DF(x0) (DEFAULT: kernel(DF(x0))
% min_res_N      double, Newton's residual ( DEFAULT: 10^-14 )
% bool_Hopf      bool, if 1 appended scalar equation during continuation
%                ||sol||=1 (DEFAULT: 0)
% 
% OUTPUT
% s              string, path of saved solutions
% x_n            Xi_vector, last validated solution
% a_lot_of_stuff STILL TO BE DETERMINED

global talkative
global use_intlab

temp_use_intval = use_intlab;
use_intlab =0;

% test of inputs
if ~isa(x0,'Xi_vector')
    error('Invalid input solution')
end
if ~isa(F_not_square,'full_problem')
    error('Invalid input system')
end
if nargin<3 || isempty(n_iter)
    n_iter = 10;
elseif ~isa(n_iter,'double')
    warning('Invalid number of iterations, reset to the default value')
    n_iter = 10;
end
if nargin<4 || isempty(h)
    h = 10^-5;
elseif ~isa(h,'double')
    warning('Invalid number of iterations, reset to the default value')
    h = 10^-5;
end
if nargin<6 || isempty(s)
    s = sprintf('simulation_of_the_%s_at_%s',datestr(now,'ddmmmyy'),datestr(now,'HHMM'));
end
if nargin<7|| isempty(min_res_N)
    min_res_N = 10^-14;
end
if nargin<8 || isempty(bool_Hopf)
    bool_Hopf =0;
end

if bool_Hopf
    if x0.size_scalar~=x0.size_vector+3
        error('If you call for a Hopf validation, the number of scalar variables should be the number of vector variables + 3')
    end
end


previous_iter.Y = [];
previous_iter.Z0 = [];
previous_iter.Z1 = [];
previous_iter.Z2 = [];

size_scalar = x0.size_scalar;
size_vector = x0.size_vector;
nodes = x0.nodes;

if F_not_square.scalar_equations.num_equations ~= size_scalar-1 &&( bool_Hopf && F_not_square.scalar_equations.num_equations ~= size_scalar-2)
    error('Invalid input problem: the number of scalar equations is inconsistent')
end

if F_not_square.vector_field.n_equations ~= size_vector
    error('Invalid input problem: the number of vector equations is inconsistent')
end

size_x = size_scalar+size_vector*(nodes*2+1);

if ~(nargin<5 || isempty(x_dot_0)) && (length(x_dot_0)~=size_x || min(size(x_dot_0)) ~= 1)
    warning('Invalid direction, default direction computed')
    x_dot_0 = [];
end


% addition of some equations for Hopf
F_old = F_not_square;
if bool_Hopf
    DF0_0 = derivative_to_matrix(derivative(F_old,x0,0));
    DF0_0(1,:)=0;
    p0 = conj(kernel(DF0_0));%[zeros(1,5),reshape(repmat(real(sum(x0.vector,2)),1,2*nodes+1)',1,[])];
    const = -Xi_vec2vec(x0).'*p0;
    F_old.scalar_equations = change_lin_coef_vector(F_old.scalar_equations,[p0;const],1);
    
    DF0_1 = derivative_to_matrix(derivative(F_old,x0,0));
    p1 = conj(kernel(DF0_1));%Xi_vec2vec(x0);
    const = -Xi_vec2vec(x0).'*p1;
    F_old.scalar_equations = change_lin_coef_vector(F_old.scalar_equations,[p1;const],2);
    F_not_square = F_old;
end

% last inout test
if nargin<5 || isempty(x_dot_0)
    DF0 = derivative_to_matrix(derivative(F_old,x0,0));
    x_dot_0 = kernel(DF0);
    if min(size(x_dot_0)) ~= 1
        error('At the computed point, the derivative is singular')
    end
end
% finished test of inputs
% construction of old system
x_dot_0 = Xi_vec2vec(symmetrise(vec2Xi_vec(x_dot_0,x0)));
const = -Xi_vec2vec(x0).'*x_dot_0;
F_old.scalar_equations = change_lin_coef_vector(F_old.scalar_equations,[x_dot_0;const],F_old.scalar_equations.number_equations_lin+1);

DF0 = derivative_to_matrix(derivative(F_old,x0,0));
Aold = inv(DF0);

% preallocate storage
step_size = zeros(n_iter,1);
norm_x = zeros(x0.size_scalar+x0.size_vector,n_iter);
Interval = zeros(2,n_iter);
Z0_iter = zeros(x0.size_scalar+x0.size_vector,n_iter);
Z1_iter = zeros(x0.size_scalar+x0.size_vector,n_iter);
Z2_iter = zeros(x0.size_scalar+x0.size_vector,n_iter);
Y_iter = zeros(x0.size_scalar+x0.size_vector,n_iter);

% test of residual in old system
[Yvector]=Y_bound_new(Aold,x0,F_old);
tries = 10;
ii = 1;
while ii<tries && max(Yvector)>10^-4
    ii = ii+1;
    x0 = Newton_2(x0,F_old,[],min_res_N);
    
    % set up some requesed elements for the validation
    DF0 = derivative_to_matrix(derivative(F_old,x0,0));
    Aold = inv(DF0);
    
    [Yvector]=Y_bound_new(Aold,x0,F_old);
end

if talkative
    fprintf('Entering validation')
end

for i =1 : n_iter
    % compute new solution with arch-length parametrisation
    x_tilde_1 = x0 + h * vec2Xi_vec(x_dot_0,x0);
    const = -Xi_vec2vec(x_tilde_1).'*x_dot_0;
    F_new = F_old;
    if bool_Hopf
        %DF0_1 = derivative_to_matrix(derivative(F_old,x0,0));
        %p1 = conj(kernel([x_dot_0.';DF0_1]));%p1 = Xi_vec2vec(x_tilde_1);
        const_p0 = -Xi_vec2vec(x_tilde_1).'*p0;
        F_new.scalar_equations = change_lin_coef_vector(F_old.scalar_equations,[p0;const_p0],1);
        const_p1 = -Xi_vec2vec(x_tilde_1).'*p1;
        F_new.scalar_equations = change_lin_coef_vector(F_new.scalar_equations,[p1;const_p1],2);
    end
    F_new.scalar_equations = change_lin_coef_vector(F_new.scalar_equations,[x_dot_0;const],F_new.scalar_equations.number_equations_lin);
    
    x1 = Newton_2(x_tilde_1,F_new,[],min_res_N);
    
    % set up some requesed elements for the validation
    DF1 = derivative_to_matrix(derivative(F_new,x1,0));
    Anew = inv(DF1);
    
    [Yvector]=Y_bound_new(Anew,x1,F_new);
    tries = 10;
    ii = 1;
    while ii<tries && max(Yvector)>10^-4
        x1 = Newton_2(x_tilde_1,F_new,[],min_res_N);
        
        % set up some requesed elements for the validation
        DF1 = derivative_to_matrix(derivative(F_new,x1,0));
        Anew = inv(DF1);
        
        [Yvector]=Y_bound_new(Anew,x1,F_new);
    end
    % validate
    [flag,Imin,Imax,previous_iter,Yvector,Z0vector,Z1vector,Z2vector,new_step] = radii_polynomials_cont_new(x0,x1,DF0,DF1,...
        F_old,F_new,previous_iter,Aold,Anew);
    
    if flag <1
        tries =1;
        while flag<1 && tries<5
            h = h* new_step;
            if talkative >0
                fprintf('could not validate the %i-th interval, h decreased to %e\n',i, h)
            end
            x_tilde_1 = x0 + h * vec2Xi_vec(x_dot_0,x0);
            const = -Xi_vec2vec(x_tilde_1).'*x_dot_0;
            %F_new = F_not_square;
            %F_new.scalar_equations = change_lin_coef_vector(F_not_square.scalar_equations,[x_dot_0;const],size_scalar);
            F_new = F_old;
            if bool_Hopf
                %DF0_1 = derivative_to_matrix(derivative(F_old,x0,0));
                %p1 = conj(kernel([x_dot_0.';DF0_1]));%p1 = Xi_vec2vec(x_tilde_1);
                const_p0 = -Xi_vec2vec(x_tilde_1).'*p0;
                F_new.scalar_equations = change_lin_coef_vector(F_old.scalar_equations,[p0;const_p0],1);
                const_p1 = -Xi_vec2vec(x_tilde_1).'*p1;
                F_new.scalar_equations = change_lin_coef_vector(F_new.scalar_equations,[p1;const_p1],2);
            end
            F_new.scalar_equations = change_lin_coef_vector(F_new.scalar_equations,[x_dot_0;const],F_new.scalar_equations.number_equations_lin);
            
            x1 = Newton_2(x_tilde_1,F_new,[],min_res_N);
            
            % set up some requesed elements for the validation
            DF1 = derivative_to_matrix(derivative(F_new,x1,0));
            Anew = inv(DF1);
            
            use_intlab = temp_use_intval;
            % validate
            [flag,Imin,Imax,previous_iter,Yvector,Z0vector,Z1vector,Z2vector,new_step] = radii_polynomials_cont_new(x0,x1,DF0,DF1,...
                F_old,F_new,previous_iter,Aold,Anew);
            
            use_intlab = 0;
        end
        if flag<1
            error('Could not validate the interval')
        end
    elseif talkative>0
        fprintf('\n ---- SUCCESS ----\n ')
        if talkative>0 && mod(i,5) ==0
            fprintf('The %i-th step has been validated successfully\n',i);
            if talkative>1
                fprintf('The stepsize is %e',h)
            end
        end
        fprintf('\n')
    end
    
    % update
    x0 = x1;
    Aold = Anew;
    F_old = F_new;
    DF0 = derivative_to_matrix(derivative(F_not_square,x0,0));
    x_dot_0 = kernel(DF0);
    x_dot_0 = Xi_vec2vec(symmetrise(vec2Xi_vec(x_dot_0,x0)));
    h = h* new_step;
    DF0 = DF1;
    
    %storage
    step_size(i+1) = h;              %zeros(n_iter,1);
    norm_x(:,i)    = vert(norm(x0)); %zeros(x0.size_scalar+x0.size_vector,n_iter);
    Interval(:,i)  = [Imin,Imax]';   %zeros(2,n_iter);
    Z0_iter(:,i)   = vert(Z0vector); %   zeros(x0.size_scalar+x0.size_vector,n_iter);
    Z1_iter(:,i)   = vert(Z1vector); %zeros(x0.size_scalar+x0.size_vector,n_iter);
    Z2_iter(:,i)   = vert(Z2vector); %zeros(x0.size_scalar+x0.size_vector,n_iter);
    Y_iter(:,i)    = vert(Yvector);  %zeros(x0.size_scalar+x0.size_vector,n_iter);
end

x_n = x1;
clear x1;
clear Aold;clear Anew; clear DF0; clear DF1;
clear h;

% save to location
save(s);

return