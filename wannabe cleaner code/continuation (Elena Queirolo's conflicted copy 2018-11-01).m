function [s, x_n] = continuation ( x0, F_not_square, n_iter, h, x_dot_0,s, min_res_N, bool_Hopf,bool_fancy_scalar, bool_saddle)
% function [s, x_n] = continuation ( x0, F_not_square, n_iter, h, x_dot_0,s,min_res_N, bool_Hopf,bool_fancy_scalar, bool_saddle)
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
%                ||sol||=1 and y_1(0)=0 (DEFAULT: 0)
% bool_fancy_scalar   bool, requesting the FIRST linear scalar equation to
%                be \int_[0,2pi] x(t) dot x'(t) dt = 0 and updated at every
%                iteration (DEFAULT: 0)
% bool_saddle    bool, if 1 check for saddle nodes after each segment
%                (numerically) and if found one, validate it (DEFAULT: 0)
%
% OUTPUT
% s              string, path of saved solutions
% x_n            Xi_vector, last validated solution
%

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
if nargin<9 || isempty(bool_fancy_scalar)
    bool_fancy_scalar = 0;
end
if nargin<10 || isempty(bool_saddle)
    bool_saddle = 0;
end

% data coming from the previous iteration
previous_iter.Y = [];
previous_iter.Z0 = [];
previous_iter.Z1 = [];
previous_iter.Z2 = [];

% some useful variables
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
    % first scalar equation: imposing x1(0) = 0
    % THIS IS THE BAD ONE
    %     p0 = zeros(size_x,1);
    %     p0(1:5) = 0;
    %     p0(6:(5+2*x0.nodes+1)) = 1;
    %     p0((5+2*x0.nodes+2):end) = 0;
    %     const_p0 = 0;
    %     F_old.scalar_equations = change_lin_coef_vector(F_old.scalar_equations,[p0;const_p0],1);
    
    % second continuation condition
    % int x .* tilde x = 1
    %     p1 = Xi_vec2vec(x0)';
    %     p1(1:5) = 0;
    %     const_p1 = -1;
    %     F_old.scalar_equations = change_lin_coef_vector(F_old.scalar_equations,[p1.';const_p1],2);
    %     F_not_square = F_old;
end

% addition of some equations for fancy_scalar_equation
if bool_fancy_scalar
    F_old.scalar_equations = fancy_scalar_condition(x0,F_old.scalar_equations);
end

% last input test
if nargin<5 || isempty(x_dot_0) % definition of x_dot if not defined yet
    DF0 = derivative_to_matrix(derivative(F_old,x0,0));
    x_dot_0 = kernel(DF0);
    if min(size(x_dot_0)) ~= 1
        error('At the computed point, the derivative is singular')
    end
end

% addition of the 10th July 2018
% if h<0
%     x_dot_0=-x_dot_0;
%     h=-h;
% end
% finished test of inputs

% construction of old system
[~,index] = max(abs(real(x_dot_0(1:x0.size_scalar))));
angle = atan( imag(x_dot_0(index))/real(x_dot_0(index)));
x_dot_0 = exp( - 1i * angle) * x_dot_0; % bringing x_dot_0 to be symmetric (by multiplication with the appropriate complex rotation)
x_dot_0 = Xi_vec2vec(symmetrise(vec2Xi_vec(x_dot_0,x0)));
const = -Xi_vec2vec(x0).'*x_dot_0;
F_old.scalar_equations = change_lin_coef_vector(F_old.scalar_equations,[x_dot_0;const],F_old.scalar_equations.number_equations_lin+1); % addition of the continuation equation

DF0 = derivative_to_matrix(derivative(F_old,x0,0));
Aold = inv(DF0);

% preallocate storage
step_size = zeros(n_iter,1);
norm_x = zeros(x0.size_scalar+x0.size_vector,n_iter+1);
Interval = zeros(2,n_iter);
Z0_iter = zeros(x0.size_scalar+x0.size_vector,n_iter);
Z1_iter = zeros(x0.size_scalar+x0.size_vector,n_iter);
Z2_iter = zeros(x0.size_scalar+x0.size_vector,n_iter);
Y_iter = zeros(x0.size_scalar+x0.size_vector,n_iter);

saddle_x0_stored =cell(1);
saddle_x1_stored =cell(1);
saddle_y0 =[];
saddle_y1 =[];
saddle_z0 =[];
saddle_z1 =[];
saddle_index_var =[];

% validating old system % need a good starting point
x0 = Newton_2(x0,F_old,[],min_res_N);
try
    if talkative >1
        fprintf('Validation of starting point\n');
    end
    radii_polynomials(x0,F_old,DF0,Aold);
catch
    if talkative >1
        fprintf('Starting point require extra precision, running Newton\n');
    end
    x0 = Newton_2(x0,F_old,[],min_res_N);
    
    % set up some requesed elements for the validation
    DF0 = derivative_to_matrix(derivative(F_old,x0,0));
    Aold = inv(DF0);
    try
        % pointwise validation
        radii_polynomials(x0,F_old,DF0,Aold);
    catch
        error('Initial point not good enough, even after Newton. Consider changing the number of nodes and or nu');
    end
end

norm_x(:,1)    = vert(norm(x0));
norm_x(1:x0.size_scalar,1) = x0.scalar;

if talkative
    fprintf('Start of continuated validation\n')
end

for i =1 : n_iter
    % compute new step (and system to solve at the new step) with a
    % predictor-corrector algorithm
    [x1,x_dot_1,F_new] =update_new(x0,h,x_dot_0,...
        F_not_square,[],min_res_N, bool_Hopf,bool_fancy_scalar);
    
    % set up some requesed elements for the validation
    DF1 = derivative_to_matrix(derivative(F_new,x1,0));
    Anew = inv(DF1);
    
    % validate
    [flag,Imin,Imax,previous_iter,Yvector,Z0vector,Z1vector,Z2vector,new_step] = radii_polynomials_cont_new(x0,x1,DF0,DF1,...
        F_old,F_new,previous_iter,Aold,Anew);
    
    if flag <1 % if the validation did not work, redo the previous steps
        tries =1;
        while flag<1 && tries<5
            h = h* new_step;
            tries = tries+1;
            if talkative >0
                fprintf('could not validate the %i-th interval, h decreased to %e\n',i, h)
            end
            % compute new step (and system to solve at the new step) with a
            % predictor-corrector algorithm
            use_intlab = 0;
            [x1,x_dot_1,F_new] =update_new(x0,h,x_dot_0,...
                F_not_square,[],min_res_N, bool_Hopf, bool_fancy_scalar);
            
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
    end
    % vadlidation succeeded or program crashed
    
    
    if bool_saddle
        
        % % % DEBUGGING
        % [numerical_check, x_prime0_int, x_prime1_int,x_prime_prime0_int,...
        % x_prime_prime1_int] = if_saddle_numerical(F_old,F_new, x0,x1,2*bool_Hopf);
        % %
        % lambda0(i,:) = [x0.scalar(2), x_prime0_int(2),x_prime_prime0_int(2)];
        % T0(i,:) = [x0.scalar(1), x_prime0_int(1),x_prime_prime0_int(1)];
        % amplitude0(i,:) = [x0.scalar(3), x_prime0_int(3),x_prime_prime0_int(3)];
        %
        % lambda1(i,:) = [x1.scalar(2), x_prime1_int(2),x_prime_prime1_int(2)];
        % T1(i,:) = [x1.scalar(1), x_prime1_int(1),x_prime_prime1_int(1)];
        % amplitude1(i,:) = [x1.scalar(3), x_prime1_int(3),x_prime_prime1_int(3)];
        saddle_variable = [];
        if bool_Hopf
            saddle_variable = 2;
        end
        if bool_Hopf
            [numerical_check] = if_saddle_numerical(F_old,F_new, x0,x1,2*bool_Hopf);
        else
            [numerical_check] = if_saddle_numerical(F_old,F_new, x0,x1);
        end
        if numerical_check
            if bool_Hopf
                [saddle_confirmed]=validation_saddle(F_old,F_new, x0,x1,2*bool_Hopf);
            else
                [saddle_confirmed]=validation_saddle(F_old,F_new, x0,x1);
            end
        else saddle_confirmed=0;
        end
        if saddle_confirmed
            % save the initial solution, plus the derivatives of the
            % scalar parts, as well as the index of the variable that
            % has a saddle node
            disp('THE SADDLE HAS BEEN CONFIRMED!')
            pause(3)
            % by the way, now we are confident of its uniqueness
            if saddle_confirmed>1
                warning('There are %i saddles in this segment', saddle_confirmed)
            end
            % %                 n_saddle = length(saddle_index_var)+1;
            % %                 saddle_x0_stored{n_saddle} = saddle_x0;
            % %                 saddle_x1_stored{n_saddle} = saddle_x1;
            % %                 saddle_y0(n_saddle,:) = horiz(y0);
            % %                 saddle_y1(n_saddle,:) = horiz(y1);
            % %                 saddle_z0(n_saddle,:) = horiz(z0);
            % %                 saddle_z1(n_saddle,:) = horiz(z1);
            % %                 saddle_index_var(n_saddle) =saddle_confirmed;
        end
    end
    
    % test on the smoothness condition 
    flag = smoothness_condition(x0, x_dot_0, x1,x_dot_1, Imin);
    if flag == 0
        error('The validated curve is not guaranteed to be smooth')
    end
    
    % check if too many or too little nodes
    % 
    % introduce 2 counters: too many, too little
    % if was already too many, new solution smaller
    % if too many, last mode(s) of new step imposed zero
    % if too little, add mode to new step asap
    % check for "ideal" number of nodes
    first_big_mode = find(log(sum(abs(x1.vector),1))>-14,1);
    ideal_small_modes = max(5,0.25*x0.nodes);
    
    if first_big_mode <= 0.9*ideal_small_modes
        add_node = floor((-first_big_mode+ideal_small_modes)*0.5+1);
        decrease_node = 0;
    elseif first_big_mode>1.2*ideal_small_modes
        decrease_node = ceil((-first_big_mode+ideal_small_modes)*0.5+1);
        add_node = 0;
    else
        add_node =0;
        decrease_node = 0;
    end
    
    if add_node
        if talkative >1
            fprintf('At iteration %i, the number of nodes is increased by %i. New number of modes %i\n', i, add_node, x0.nodes+add_node);
        end
        x1 = reshape_Xi(x1, x1.nodes+add_node);
        x_dot_1 = Xi_vec2vec(reshape_Xi(vec2Xi_vec(x_dot_1,x0), ...
            x0.nodes+add_node));
        F_new.scalar_equations = reshape(F_new.scalar_equations, x0.nodes+add_node);
        F_not_square.scalar_equations = reshape(F_not_square.scalar_equations, x0.nodes+add_node);
        DF1 = derivative_to_matrix(derivative(F_new,x1,0));
        Anew = inv(DF1);
    end
    if decrease_node 
        if talkative >1
            fprintf('At iteration %i, the number of nodes is decreased by %i. New number of modes %i\n', i, abs(decrease_node), x0.nodes+decrease_node);
        end
        % necessary temporary elements for "one point validation"
        x1_big = x1;
        Anew_big = Anew;
        DF1_big = DF1;
        F_new_big = F_new;
        x1 = reshape_Xi(x1, x1.nodes+decrease_node);
        x_dot_1 = Xi_vec2vec(reshape_Xi(vec2Xi_vec(x_dot_1,x0),...
            x0.nodes+decrease_node));
        F_new.scalar_equations = reshape(F_new.scalar_equations, x0.nodes+decrease_node);
        F_not_square.scalar_equations = reshape(F_not_square.scalar_equations, x0.nodes+decrease_node);
        % following node with less modes
        
        % one point validation 
        x1_small = reshape_Xi(x1, x0.nodes);
        F_new_small = F_new;
        F_not_square_small = F_not_square ;
        F_new_small.scalar_equations = reshape(F_new.scalar_equations, x0.nodes);
        F_not_square_small.scalar_equations = reshape(F_not_square.scalar_equations, x0.nodes);
        DF1_small = derivative_to_matrix(derivative(F_new_small,x1_small,0));
        Anew_small = inv(DF1_small);
        
        use_intlab = temp_use_intval;
        % validate
        [flag] = radii_polynomials_cont_new(x1_big,x1_small,DF1_big,DF1_small,...
            F_new_big,F_new_small,[],Anew_big,Anew_small);
        
        use_intlab = 0;
        if flag==0
            error('When decreasing number of modes something went wrong')
        end
        
        DF1 = derivative_to_matrix(derivative(F_new,x1,0));
        Anew = inv(DF1);
    end
    
    
    if talkative>0
        % some talking to the user, when requested
        fprintf('\n ---- SUCCESS ----\n ')
        if talkative>1
            fprintf('The %i-th interval of validation is [%e, %e]\n',i,Imin,Imax);
        end
        if talkative>0 && mod(i,5) ==0
            fprintf('The %i-th step has been validated successfully\n',i);
            if talkative>1
                fprintf('The stepsize is %e\n',h)
            end
        end
        fprintf('\n')
    end
    
    % iteration
    x_dot_0 = x_dot_1;
    Aold = Anew;
    F_old = F_new;
    h = h* new_step;
    DF0 = DF1;
    x0 = x1;
    
    %storage
    step_size(i+1) = h;
    norm_x(:,i+1)    = vert(norm(x0));
    norm_x(1:x0.size_scalar,i+1) = x0.scalar;
    Interval(:,i)  = [Imin,Imax]';
    Z0_iter(:,i)   = vert(Z0vector);
    Z1_iter(:,i)   = vert(Z1vector);
    Z2_iter(:,i)   = vert(Z2vector);
    Y_iter(:,i)    = vert(Yvector);
end

x_n = x1;
x_dot_n = x_dot_1;

% save to location the elements that are interesting
if bool_saddle && ~isempty(saddle_index_var)
    save(s, 'x_n', 'x_dot_n', 'step_size', 'norm_x', 'Interval', 'Y_iter', 'Z0_iter', 'Z1_iter', 'Z2_iter', ...
        'saddle_x0_stored','saddle_x1_stored','saddle_y0','saddle_y1','saddle_z0','saddle_z','saddle_index_var');
else
    save(s, 'x_n', 'x_dot_n','step_size', 'norm_x', 'Interval', 'Y_iter', 'Z0_iter', 'Z1_iter', 'Z2_iter');
end

return