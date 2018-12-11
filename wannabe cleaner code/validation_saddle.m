function [saddle_confirmed]=validation_saddle(F_old,F_new, x0,x1,num_variable)
global use_intlab 
global talkative
global refinement_saddle
global rescaling_saddle 
temp_use_intlab = use_intlab;
use_intlab = 0;

if ~exist('num_variable','var')|| ~exist('num_variable','var')
    num_variable=1:x0.size_scalar;
end

[numerical_check,~,~,~,~,num_variable] = if_saddle_numerical(F_old,F_new, x0,x1, num_variable);
 
if ~numerical_check
    saddle_confirmed = 0;
    return
end
 
saddle_confirmed = 0;

if isempty(refinement_saddle)
    n_intervals= 300;
else
    n_intervals = refinement_saddle;
end
 
% 
% [numerical_check, x_prime0, x_prime1,x_prime_prime0,x_prime_prime1] = if_saddle_numerical(F_old,F_new, x0,x1);
% step = (1/n_intervals)*(x1 - x0);
% x_old_step = x0;
% lin_coef_step{1} = (1/n_intervals)*(F_new.scalar_equations.linear_coef{1} - F_old.scalar_equations.linear_coef{1});
% lin_coef_step{2} = (1/n_intervals)*(F_new.scalar_equations.linear_coef{2} - F_old.scalar_equations.linear_coef{2});
% lin_coef_step{3} = (1/n_intervals)*(F_new.scalar_equations.linear_coef{3} - F_old.scalar_equations.linear_coef{3});
% F_old_step = F_old;
% numerical_check_step = zeros(n_intervals,1);
% y0_step = zeros(n_intervals,2);
% z0_step = zeros(n_intervals,2);
% y1_step = zeros(n_intervals,2);
% z1_step = zeros(n_intervals,2);
% for i = 1: n_intervals
%     x_new_step = x_old_step + step;
%     F_new_step = F_old_step;
%     
%     F_new_step.scalar_equations.linear_coef{1} = (F_old_step.scalar_equations.linear_coef{1} + lin_coef_step{1});
%     F_new_step.scalar_equations.linear_coef{2} = (F_old_step.scalar_equations.linear_coef{2} + lin_coef_step{2});
%     F_new_step.scalar_equations.linear_coef{3} = (F_old_step.scalar_equations.linear_coef{3} + lin_coef_step{3});
%     
%     [numerical_check_step(i),x_prime0t, x_prime1t,x_prime_prime0t,x_prime_prime1t] = if_saddle_numerical(F_old_step,F_new_step, x_old_step,x_new_step);
%     y0_step(i,:) = x_prime0t(1:2);
%     y1_step(i,:) = x_prime1t(1:2);
%     z0_step(i,:) = x_prime_prime0t(1:2);
%     z1_step(i,:) = x_prime_prime1t(1:2);
%     if numerical_check_step(i)
%         %[~, x_prime0, x_prime1,x_prime_prime0,x_prime_prime1] = if_saddle_numerical(F_old,F_new, x0,x1);
%         [saddle_x0, saddle_x1,y0_point,y1_point,z0_point,z1_point] = precise_computation(F_old_step, F_new_step, x_old_step,x_new_step, x_prime0, x_prime1,x_prime_prime0,x_prime_prime1);
%         y0_step(i,:) = y0_point(1:2);
%         y1_step(i,:) = y1_point(1:2);
%         z0_step(i,:) = z0_point(1:2);
%         z1_step(i,:) = z1_point(1:2);
%     end
%     x_old_step = x_new_step;
%     F_old_step = F_new_step;
% end
% pause(1)
 
 
 
 
 
for i = 1: n_intervals
    [x0_int, x1_int, F0_int, F1_int] = automatic_refinement(F_old, F_new, x0, x1, n_intervals, i); % done
    [numerical_check, x_prime0_int, x_prime1_int,x_prime_prime0_int,...
        x_prime_prime1_int] = if_saddle_numerical(F0_int, F1_int, x0_int, x1_int); %done
     
    if ~numerical_check
        continue
    end
     
    [saddle_F0_int, saddle_F1_int] = bigger_saddle_system(F0_int, F1_int, x0_int, x1_int); %done
    
    % DEBUG RYCHKOV
    %num_variable = 2;
    rescale_derivative = rescaling_saddle(1);
    rescale_second_der = rescaling_saddle(end);
    saddle_F0_int_loc = rescale(saddle_F0_int, x0.size_scalar+num_variable, rescale_derivative);
    saddle_F1_int_loc = rescale(saddle_F1_int, x0.size_scalar+num_variable, rescale_derivative);
    saddle_F0_int_loc = rescale(saddle_F0_int_loc, 2*x0.size_scalar+num_variable, rescale_second_der);
    saddle_F1_int_loc = rescale(saddle_F1_int_loc, 2*x0.size_scalar+num_variable, rescale_second_der);
    x_prime0_int(num_variable) = rescale_derivative * x_prime0_int(num_variable);
    x_prime_prime0_int(num_variable) = rescale_second_der *x_prime_prime0_int(num_variable);
    x_prime1_int(num_variable) = rescale_derivative * x_prime1_int(num_variable);
    x_prime_prime1_int(num_variable) = rescale_second_der *x_prime_prime1_int(num_variable);
    
    saddle_x0_int = compute_saddle(saddle_F0_int_loc,x0_int,x_prime0_int,x_prime_prime0_int); % done
    saddle_x1_int = compute_saddle(saddle_F1_int_loc,x1_int,x_prime1_int,x_prime_prime1_int); % done
     
    % %DEBUG
    % figure(1); hold on; plot(i,saddle_x0_int.scalar(x0.size_scalar+ num_variable),'*')
    % figure(2); hold on; plot(i,saddle_x0_int.scalar(2*x0.size_scalar+ num_variable),'*')
    % continue
    % %debug
     
    [bool_saddle, z_crosses]=check_if_saddle_possible(saddle_x0_int,saddle_x1_int,num_variable); % done
    if bool_saddle && z_crosses
        error('what to do now?')
    end
    if bool_saddle
        use_intlab = temp_use_intlab;
        [delta, bool_saddle_confirmed] = bound_refinement(saddle_x0_int, saddle_x1_int, saddle_F0_int_loc, saddle_F1_int_loc,num_variable); % done
        use_intlab = 0;
        if bool_saddle_confirmed
        delta =3; bool_saddle_confirmed=0;
        end
        if bool_saddle_confirmed
            use_intlab = temp_use_intlab;
            saddle_confirmed = saddle_validation(saddle_x0_int, saddle_x1_int, saddle_F0_int_loc, saddle_F1_int_loc,num_variable);
            use_intlab = temp_use_intlab;
            if saddle_confirmed==0
                error('Adding intlab is problematic')
            end
            return
        end
        n_intervals_ref = ceil(1/delta)+1; %ceil(norm(norm(saddle_x0_int-saddle_x1_int))/delta)+1;
        for j = 1: n_intervals_ref
            [saddle_x0_small_int, saddle_x1_small_int, F0_small_int,...
                F1_small_int] = automatic_refinement(saddle_F0_int_loc, ...
                saddle_F1_int_loc, saddle_x0_int, saddle_x1_int, n_intervals_ref, j);% done
            saddle_x0_small_int = Newton_2(saddle_x0_small_int,F0_small_int);
            saddle_x1_small_int = Newton_2(saddle_x1_small_int,F1_small_int);
            bool_saddle_small_int=check_if_saddle_possible(saddle_x0_small_int,saddle_x1_small_int,num_variable); % done
            if bool_saddle_small_int
                use_intlab = temp_use_intlab;
                bool_saddle_tot = saddle_validation(saddle_x0_small_int, saddle_x1_small_int, F0_small_int, F1_small_int,num_variable);
                use_intlab = 0;
                if bool_saddle_tot
                    saddle_confirmed = saddle_confirmed+1;
                    use_intlab = temp_use_intlab;
                    return
                else
                    warning('Failure to validate')
                end
            end
        end
    else
        if talkative >2
            disp('False positive, not a saddle');
        end
    end
     
end
use_intlab = temp_use_intlab;
end
 
% functions to do
% [x0_int, x1_int, F0_int, F1_int] = automatic_refinement(F_old, F_new, x0,
% x1, n_intervals, i); % done
% saddle_x0_int = compute_saddle(saddle_F0_int,x0_int); % done
% [bool_saddle, z_crosses]=check_if_saddle_possible(x0_int,x1_int); % done
% [delta, bool_saddle_confirmed] = bound_refinement(saddle_x0_int,
% saddle_x1_int, saddle_F0_int, saddle_F1_int); % done
% bool_saddle_tot = saddle_validation(saddle_x0_small_int, saddle_x1_small_int, F0_small_int, F1_small_int);
 
 
function [x0_int, x1_int, F0_int, F1_int] = automatic_refinement(F_old, F_new, x0, x1, n_intervals, i)
% function [x0_int, x1_int, F0_int, F1_int] = automatic_refinement(F_old, F_new, x0, x1, n_intervals, i)
% 
% split an interval into i subintervals and return the numerical solutions
% and the vector fields of the subintervals
 
x_step = (1/n_intervals)*(x1 - x0);
lin_coef_step{1} = (1/n_intervals)*(F_new.scalar_equations.linear_coef{1} - F_old.scalar_equations.linear_coef{1});
lin_coef_step{2} = (1/n_intervals)*(F_new.scalar_equations.linear_coef{2} - F_old.scalar_equations.linear_coef{2});
lin_coef_step{3} = (1/n_intervals)*(F_new.scalar_equations.linear_coef{3} - F_old.scalar_equations.linear_coef{3});
 
x0_int = x0 + (i-1) * x_step;
x1_int = x0 + i * x_step;
 
F0_int = F_old;
F1_int = F_old;
 
F0_int.scalar_equations.linear_coef{1} = F_old.scalar_equations.linear_coef{1} + (i-1) *lin_coef_step{1};
F0_int.scalar_equations.linear_coef{2} = F_old.scalar_equations.linear_coef{2} + (i-1) *lin_coef_step{2};
F0_int.scalar_equations.linear_coef{3} = F_old.scalar_equations.linear_coef{3} + (i-1) *lin_coef_step{3};
 
F1_int.scalar_equations.linear_coef{1} = F_old.scalar_equations.linear_coef{1} + i *lin_coef_step{1};
F1_int.scalar_equations.linear_coef{2} = F_old.scalar_equations.linear_coef{2} + i *lin_coef_step{2};
F1_int.scalar_equations.linear_coef{3} = F_old.scalar_equations.linear_coef{3} + i *lin_coef_step{3};
 
end
 
function saddle_x = compute_saddle(saddle_F,x,y,z)
% function saddle_x = compute_saddle(saddle_F,x)
%
% inputs: big system and small solution
% output: big solution
y = vec2Xi_vec(y,x);
z = vec2Xi_vec(z,x);
saddle_x0 = Xi_vector([x.scalar,y.scalar,z.scalar],[x.vector;y.vector;z.vector]);
saddle_x = Newton_2(saddle_x0, saddle_F);
end
 
function [bool_saddle, z_crosses]=check_if_saddle_possible(x0,x1,num_variable)
% function [bool_saddle, z_crosses]=check_if_saddle_possible(x0_int,x1_int)
%
% inputs: 2 big solutions
% output: 2 bools, check if the first derivative is crossing, and if the
% second is crossing too the zero axes
size_scalar = x0.size_scalar /3;
if size_scalar ~=ceil(size_scalar)
    error('inputs off')
end
 
y0 = (x0.scalar(size_scalar+1:size_scalar*2));
y0 = y0(num_variable);
y1 = (x1.scalar(size_scalar+1:size_scalar*2));
y1 = y1(num_variable);

z0 = (x0.scalar(size_scalar*2+1:size_scalar*3));
z0 = z0(num_variable);
z1 = (x1.scalar(size_scalar*2+1:size_scalar*3));
z1 = z1(num_variable);

index_saddle = find( y0.*y1<0);
if ~isempty(index_saddle)
    bool_saddle=1;
else
    bool_saddle=0;
    z_crosses =0;
    return
end
if any(z0(index_saddle).*z1(index_saddle)<0)
    z_crosses=1;
else
    z_crosses =0;
end
 
end
 
 
function [delta, bool_saddle_confirmed] = bound_refinement(saddle_x0, saddle_x1, saddle_F0, saddle_F1,num_variable)
% validate the interval and define the maximum length of the interval to be
% able to confirm that z doesn't change sign in that interval
bool_saddle_confirmed =0;
 
DF0_saddle = derivative_to_matrix(derivative(saddle_F0,saddle_x0,0));
Aold_saddle = inv(DF0_saddle);
DF1_saddle = derivative_to_matrix(derivative(saddle_F1,saddle_x1,0));
Anew_saddle = inv(DF1_saddle);
[flag,Imin,~,~,~,~,~,~,~,Ycont] = radii_polynomials_cont_new(saddle_x0,saddle_x1,DF0_saddle,DF1_saddle,...
    saddle_F0,saddle_F1,[],Aold_saddle,Anew_saddle);
 
size_scalar = saddle_x0.size_scalar/3;
z0 = (saddle_x0.scalar(size_scalar*2+1:size_scalar*3));
z1 = (saddle_x1.scalar(size_scalar*2+1:size_scalar*3));
z_bound = min(min(abs(z0(num_variable))),min(abs(z1(num_variable))));
if flag>0 && Imin < z_bound
    % validated already
    bool_saddle_confirmed =1 ;
    delta = 1;
    return
    %error('It happened!') % still to code, laziness is a bad habit
end
Y_extrema_max = max(Ycont(:,4));
Y_s_max = max(Ycont(:,2));
 
if z_bound>Y_extrema_max%%% NEW LINE
delta = sqrt( (z_bound - Y_extrema_max)/Y_s_max);
else%%% NEW LINE
    delta = 1;%%% NEW LINE
end%%% NEW LINE
 
end
 
 
function saddle_confirmed = saddle_validation(saddle_x0, saddle_x1, saddle_problem0, saddle_problem1,num_variable)
 
global talkative
 
DF0_saddle = derivative_to_matrix(derivative(saddle_problem0,saddle_x0,0));
Aold_saddle = inv(DF0_saddle);
DF1_saddle = derivative_to_matrix(derivative(saddle_problem1,saddle_x1,0));
Anew_saddle = inv(DF1_saddle);
[flag,Imin_sad] = radii_polynomials_cont_new(saddle_x0,saddle_x1,DF0_saddle,DF1_saddle,...
    saddle_problem0,saddle_problem1,[],Aold_saddle,Anew_saddle);
if flag <=0
    error('Could not validate saddle')
end
size_scalar = saddle_x0.size_scalar/3;
y0_point = (saddle_x0.scalar(size_scalar+1:size_scalar*2));
y1_point = (saddle_x1.scalar(size_scalar+1:size_scalar*2));
 
z0_point = (saddle_x0.scalar(size_scalar*2+1:size_scalar*3));
z1_point = (saddle_x1.scalar(size_scalar*2+1:size_scalar*3));
 
% check if validation confirmed saddle
% if any(abs(z0_point(num_variable))<Imin_sad) || ...
%         any(abs(z1_point(num_variable))<Imin_sad) ||...
%     any(abs(y0_point(num_variable))<Imin_sad) || ...
%         any(abs(y1_point(num_variable))<Imin_sad)
%     error('Bounds still too big!')
% end
if any(abs(z0_point(num_variable))<Imin_sad) || ...
        any(abs(z1_point(num_variable))<Imin_sad) %||...
    %any(abs(y0_point(num_variable))<Imin_sad) || ...
    %    any(abs(y1_point(num_variable))<Imin_sad)
    error('Bounds for z still too big! Solution by rescaling')
end
Imin_sad_y0 = Imin_sad;
Imin_sad_y1 = Imin_sad;
if any(abs(y0_point(num_variable))<Imin_sad_y0) 
    [~,Imin_sad_y0]=radii_polynomials(saddle_x0,saddle_problem0,DF0_saddle, Aold_saddle);
end
if any(abs(y1_point(num_variable))<Imin_sad_y1)
    [~,Imin_sad_y1]=radii_polynomials(saddle_x1,saddle_problem1,DF1_saddle, Anew_saddle);
end
if any(abs(y0_point(num_variable))<Imin_sad_y0) ||...
    any(abs(y1_point(num_variable))<Imin_sad_y1)
    error('Bounds dor y still too big! Consider rescaling?')
end

y0 = midrad(y0_point,Imin_sad);
y1 = midrad(y1_point,Imin_sad);
 
z0 = midrad(z0_point,Imin_sad);
z1 = midrad(z1_point,Imin_sad);
 
saddle_confirmed = 0;
index_saddle = find(y1(num_variable).*y0(num_variable)<0);
if isempty(index_saddle)
    [flag_old,Imin_old]=radii_polynomials(saddle_x0,saddle_problem0,DF0_saddle,Aold_saddle);
    [flag_new,Imin_new]=radii_polynomials(saddle_x1,saddle_problem1,DF1_saddle,Anew_saddle);
    if flag_old*flag_new ==0
        error('NOOOOOOO');
    end
    I_min = min(Imin_old,Imin_new);
     
    y0 = midrad(y0_point,I_min);
    y1 = midrad(y1_point,I_min);
     
    index_saddle = find(y1.*y0<0);
    if isempty(index_saddle)
        if isempty(find(y1.*y0>0, 1))
            error('nope')
        else
            return
        end
    end
end
index_saddle= find(z0(num_variable(index_saddle)).*z1(num_variable(index_saddle))>0);
if length(index_saddle)==1
    saddle_confirmed = 1;
    if talkative>0
        fprintf('\n       A saddle node was confirmed and is getting stored\n')
    end
elseif  isempty(index_saddle)
    saddle_confirmed = 0;
    if talkative>0
        fprintf('\n       A saddle node was not confirmed and will not be stored :(\n')
    end
else
    warning('Something is weird here: multiple variables that have a saddle node')
end
end
 
%%%%%%%%%%%% NOT USED ANYMORE
function [saddle_confirmed]=validation_saddle_nope(F_old,F_new, x0,x1)
% function [saddle_confirmed]=validation_saddle(F_old,F_new, x0,x1)
%
% validation of a saddle node between x0 and x1
% outline:
% numerical check
% construction of bigger system
% Newton
% validation of bigger system
% end checks
% storage
%
% INPUT
% F_old, F_new      full_problems, square problems, including continuation
%                   equation
% x0, x1            Xi_vectors, validated solutions to F_old and F_new
%                   respectively
% OUTPUT
% saddle_confirmed  bool, 1 if the saddle bifurcation was validated, 0
%                   otherwise
 
% other outputs to add
% saddle_index_var saddle_x0 saddle_x1 y0 y1 z0 z1
% %
% % % new pseudocode:
% % % numerical check
% % % test on [x0,x1]
% % % if radmin > abs(y) or abs(z) AND || x0  - x1 || > delta
% % %   if numerical check (x0, midpoint)
% % %    test on [x0, mid_point]
% % %   else if numerical check (midpoint, x1)
% % %    test on [mid_point, x0]
% % %   else
% % %     wtf
% % %    end
% % % else
% % %    impossible validation
% % % end
 
global talkative
 
% numerical check
[numerical_check, x_prime0, x_prime1,x_prime_prime0,x_prime_prime1] = if_saddle_numerical(F_old,F_new, x0,x1);
 
if ~numerical_check
    saddle_confirmed = 0;
    return
end
% [saddle_confirmed,bounds_too_big] = actual_validation(F_old, F_new, x0,x1, x_prime0, x_prime1,x_prime_prime0,x_prime_prime1);
 
 
% % if saddle_confirmed ==1
% %     return
% % elseif bounds_too_big ==1 && norm(norm(x0 - x1))>10^-9
% %     % divide segment in half
%     x_mid = 0.5*(x0+x1);
%     F_mid = F_old;
%     F_mid.scalar_equations.linear_coef{1} = 0.5*(F_old.scalar_equations.linear_coef{1} + F_new.scalar_equations.linear_coef{1});
%     F_mid.scalar_equations.linear_coef{2} = 0.5*(F_old.scalar_equations.linear_coef{2} + F_new.scalar_equations.linear_coef{2});
%     F_mid.scalar_equations.linear_coef{3} = 0.5*(F_old.scalar_equations.linear_coef{3} + F_new.scalar_equations.linear_coef{3});
%
%     [numerical_check1] = if_saddle_numerical(F_old,F_mid, x0,x_mid);
%     [numerical_check2] = if_saddle_numerical(F_mid,F_new, x_mid,x1);
%     if numerical_check1 ==0 && numerical_check2 ==0
%         error('something went weird');
%     elseif numerical_check1
%         [saddle_confirmed] = actual_validation(F_old, F_mid, x0,x_mid, x_prime0, x_prime1,x_prime_prime0,x_prime_prime1);
%         %[saddle_confirmed]=validation_saddle(F_old,F_mid, x0,x_mid);
%     else
%         [saddle_confirmed] = actual_validation(F_mid, F_new, x_mid,x1, x_prime0, x_prime1,x_prime_prime0,x_prime_prime1);
%         %[saddle_confirmed]=validation_saddle(F_mid,F_new, x_mid,x1);
%     end
% % end
z0 = x_prime_prime0(1:x0.size_scalar);
z1 = x_prime_prime1(1:x0.size_scalar);
bound = min( min ( abs( z0), abs(z1)))/2;
 
%maximum_length = test_validation(bound,F_old, F_new, x0,x1, x_prime0, x_prime1,x_prime_prime0,x_prime_prime1);
 
n_intervals = 30;%ceil(norm(norm(x0-x1))/maximum_length)+1;
 
step = (1/n_intervals)*(x1 - x0);
x_old_step = x0;
lin_coef_step{1} = (1/n_intervals)*(F_new.scalar_equations.linear_coef{1} - F_old.scalar_equations.linear_coef{1});
lin_coef_step{2} = (1/n_intervals)*(F_new.scalar_equations.linear_coef{2} - F_old.scalar_equations.linear_coef{2});
lin_coef_step{3} = (1/n_intervals)*(F_new.scalar_equations.linear_coef{3} - F_old.scalar_equations.linear_coef{3});
F_old_step = F_old;
numerical_check_step = zeros(n_intervals,1);
y0_step = zeros(n_intervals,2);
z0_step = zeros(n_intervals,2);
y1_step = zeros(n_intervals,2);
z1_step = zeros(n_intervals,2);
for i = 1: n_intervals
    x_new_step = x_old_step + step;
    F_new_step = F_old_step;
     
    F_new_step.scalar_equations.linear_coef{1} = (F_old_step.scalar_equations.linear_coef{1} + lin_coef_step{1});
    F_new_step.scalar_equations.linear_coef{2} = (F_old_step.scalar_equations.linear_coef{2} + lin_coef_step{2});
    F_new_step.scalar_equations.linear_coef{3} = (F_old_step.scalar_equations.linear_coef{3} + lin_coef_step{3});
     
    [numerical_check_step(i),x_prime0t, x_prime1t,x_prime_prime0t,x_prime_prime1t] = if_saddle_numerical(F_old_step,F_new_step, x_old_step,x_new_step);
    y0_step(i,:) = x_prime0t(1:2);
    y1_step(i,:) = x_prime1t(1:2);
    z0_step(i,:) = x_prime_prime0t(1:2);
    z1_step(i,:) = x_prime_prime1t(1:2);
    if numerical_check_step(i)
        %[~, x_prime0, x_prime1,x_prime_prime0,x_prime_prime1] = if_saddle_numerical(F_old,F_new, x0,x1);
        [saddle_x0, saddle_x1,y0_point,y1_point,z0_point,z1_point] = precise_computation(F_old_step, F_new_step, x_old_step,x_new_step, x_prime0, x_prime1,x_prime_prime0,x_prime_prime1);
        y0_step(i,:) = y0_point(1:2);
        y1_step(i,:) = y1_point(1:2);
        z0_step(i,:) = z0_point(1:2);
        z1_step(i,:) = z1_point(1:2);
    end
    x_old_step = x_new_step;
    F_old_step = F_new_step;
end
pause(1)
end
 
function maximum_step = test_validation(bound,F_old, F_new, x0,x1, x_prime0, x_prime1,x_prime_prime0,x_prime_prime1)
 
[saddle_problem0,saddle_problem1] = bigger_saddle_system(F_old,F_new,x0, x1);
 
y0 = vec2Xi_vec(x_prime0,x0);
y1 = vec2Xi_vec(x_prime1,x1);
z0 = vec2Xi_vec(x_prime_prime0,x0);
z1 = vec2Xi_vec(x_prime_prime1,x1);
saddle_x0 = Xi_vector([x0.scalar,y0.scalar,z0.scalar],[x0.vector;y0.vector;z0.vector]);
saddle_x1 = Xi_vector([x1.scalar,y1.scalar,z1.scalar],[x1.vector;y1.vector;z1.vector]);
 
saddle_x0=Newton_2(saddle_x0,saddle_problem0);
saddle_x1=Newton_2(saddle_x1,saddle_problem1);
 
% reset x to initial values, it cannot be changed from what it was
saddle_x0.scalar(1:x0.size_scalar) = x0.scalar;
saddle_x0.vector(1:x0.size_vector,:) = x0.vector;
saddle_x1.scalar(1:x1.size_scalar) = x1.scalar;
saddle_x1.vector(1:x1.size_vector,:) = x1.vector;
 
% VALIDATION
 
DF0_saddle = derivative_to_matrix(derivative(saddle_problem0,saddle_x0,0));
Aold_saddle = inv(DF0_saddle);
DF1_saddle = derivative_to_matrix(derivative(saddle_problem1,saddle_x1,0));
Anew_saddle = inv(DF1_saddle);
[flag,Imin,~,~,~,~,~,~,~,Ycont] = radii_polynomials_cont_new(saddle_x0,saddle_x1,DF0_saddle,DF1_saddle,...
    saddle_problem0,saddle_problem1,[],Aold_saddle,Anew_saddle);
if flag>0 && Imin < bound
    % validated already
    pause(1)
end
Y_extrema_max = max(Ycont(:,4));
Y_s_max = max(Ycont(:,2));
 
maximum_step = sqrt( (bound - Y_extrema_max)/Y_s_max);
 
end
 
function [saddle_x0, saddle_x1,y0_point,y1_point,z0_point,z1_point] = precise_computation(F_old, F_new, x0,x1, x_prime0, x_prime1,x_prime_prime0,x_prime_prime1)
% outline
% construction of bigger system
% Newton
% validation of bigger system
% end checks
 
saddle_confirmed = 0;
bounds_too_big = 0;
 
% setting for validation
% construct bigger system
[saddle_problem0,saddle_problem1, DxxHyy] = bigger_saddle_system(F_old,F_new,x0, x1); % DxxHyy still to debug: dimensions of linear coefficients not matching
 
% not in use anymore TEMPORARY
% index_old_x = 2; % index_old_x is the index of the scalar variable w.r.t. which we have the bifurcation
 
y0 = vec2Xi_vec(x_prime0,x0);
y1 = vec2Xi_vec(x_prime1,x1);
z0 = vec2Xi_vec(x_prime_prime0,x0);
z1 = vec2Xi_vec(x_prime_prime1,x1);
saddle_x0 = Xi_vector([x0.scalar,y0.scalar,z0.scalar],[x0.vector;y0.vector;z0.vector]);
saddle_x1 = Xi_vector([x1.scalar,y1.scalar,z1.scalar],[x1.vector;y1.vector;z1.vector]);
 
% solve for z variables
% by direct computation - not necessary anymore
% z0 = second_archlength_der(F_old,F_old,F_new,DxxHyy,x0,y0,saddle_x0);
% z1 = second_archlength_der(F_new,F_old,F_new,DxxHyy,x1,y1,saddle_x1);
%
% z0 = vec2Xi_vec(z0,x0);
% z1 = vec2Xi_vec(z1,x1);
%
% saddle_x0 = Xi_vector([x0.scalar,y0.scalar,z0.scalar],[x0.vector;y0.vector;z0.vector]);
% saddle_x1 = Xi_vector([x1.scalar,y1.scalar,z1.scalar],[x1.vector;y1.vector;z1.vector]);
 
saddle_x0=Newton_2(saddle_x0,saddle_problem0);
saddle_x1=Newton_2(saddle_x1,saddle_problem1);
y0_point = (saddle_x0.scalar(x0.size_scalar+1:x0.size_scalar*2));
y1_point = (saddle_x1.scalar(x0.size_scalar+1:x0.size_scalar*2));
 
z0_point = (saddle_x0.scalar(x0.size_scalar*2+1:x0.size_scalar*3));
z1_point = (saddle_x1.scalar(x0.size_scalar*2+1:x0.size_scalar*3));
end
 
 
function [saddle_confirmed, bounds_too_big] = actual_validation(F_old, F_new, x0,x1, x_prime0, x_prime1,x_prime_prime0,x_prime_prime1)
% outline
% construction of bigger system
% Newton
% validation of bigger system
% end checks
 
saddle_confirmed = 0;
bounds_too_big = 0;
 
% setting for validation
% construct bigger system
[saddle_problem0,saddle_problem1, DxxHyy] = bigger_saddle_system(F_old,F_new,x0, x1); % DxxHyy still to debug: dimensions of linear coefficients not matching
 
% not in use anymore TEMPORARY
% index_old_x = 2; % index_old_x is the index of the scalar variable w.r.t. which we have the bifurcation
 
y0 = vec2Xi_vec(x_prime0,x0);
y1 = vec2Xi_vec(x_prime1,x1);
z0 = vec2Xi_vec(x_prime_prime0,x0);
z1 = vec2Xi_vec(x_prime_prime1,x1);
saddle_x0 = Xi_vector([x0.scalar,y0.scalar,z0.scalar],[x0.vector;y0.vector;z0.vector]);
saddle_x1 = Xi_vector([x1.scalar,y1.scalar,z1.scalar],[x1.vector;y1.vector;z1.vector]);
 
% solve for z variables
% by direct computation - not necessary anymore
% z0 = second_archlength_der(F_old,F_old,F_new,DxxHyy,x0,y0,saddle_x0);
% z1 = second_archlength_der(F_new,F_old,F_new,DxxHyy,x1,y1,saddle_x1);
%
% z0 = vec2Xi_vec(z0,x0);
% z1 = vec2Xi_vec(z1,x1);
%
% saddle_x0 = Xi_vector([x0.scalar,y0.scalar,z0.scalar],[x0.vector;y0.vector;z0.vector]);
% saddle_x1 = Xi_vector([x1.scalar,y1.scalar,z1.scalar],[x1.vector;y1.vector;z1.vector]);
 
saddle_x0=Newton_2(saddle_x0,saddle_problem0);
saddle_x1=Newton_2(saddle_x1,saddle_problem1);
 
% reset x to initial values, it cannot be changed from what it was
saddle_x0.scalar(1:x0.size_scalar) = x0.scalar;
saddle_x0.vector(1:x0.size_vector,:) = x0.vector;
saddle_x1.scalar(1:x1.size_scalar) = x1.scalar;
saddle_x1.vector(1:x1.size_vector,:) = x1.vector;
 
% VALIDATION
 
% % with a particular Newton method
% iter_Newton = 0;
% while iter_Newton < 30 && max(norm(apply(saddle_problem0,saddle_x0)))>min_res_N
%     saddle_x0 = saddle_x0 - vec2Xi_vec( derivative_to_matrix(derivative(saddle_problem0,saddle_x0,0))\ Xi_vec2vec(apply(saddle_problem0,saddle_x0)) , saddle_x0);
%     % saddle_x0_vec = Xi_vec2vec(saddle_x0);
%     % saddle_x = saddle_x0_vec - derivative_to_matrix(derivative(saddle_problem0,saddle_x0,0))\ Xi_vec2vec(apply(saddle_problem0,saddle_x0));
%     % saddle_x(index_old_x) = saddle_x0_vec(index_old_x);
%     % saddle_x0 = vec2Xi_vec(saddle_x,saddle_x0);
% end
% if max(norm(apply(saddle_problem0,saddle_x0)))>min_res_N
%     error('Newton for the saddle node validation at x0 did not converge')
% end
% iter_Newton = 0;
% while iter_Newton < 30 && max(norm(apply(saddle_problem1,saddle_x1)))>min_res_N
%     saddle_x1_vec = Xi_vec2vec(saddle_x1);
%     saddle_x = saddle_x1_vec - derivative_to_matrix(derivative(saddle_problem1,saddle_x1,0))\ Xi_vec2vec(apply(saddle_problem1,saddle_x1));
%     saddle_x(index_old_x) = saddle_x1_vec(index_old_x);
%     saddle_x1 = vec2Xi_vec(saddle_x,saddle_x1);
% end
% if max(norm(apply(saddle_problem1,saddle_x1)))>min_res_N
%     error('Newton for the saddle node validation at x1 did not converge')
% end
 
global talkative
 
% validate segment
 
DF0_saddle = derivative_to_matrix(derivative(saddle_problem0,saddle_x0,0));
Aold_saddle = inv(DF0_saddle);
DF1_saddle = derivative_to_matrix(derivative(saddle_problem1,saddle_x1,0));
Anew_saddle = inv(DF1_saddle);
[flag,Imin_sad] = radii_polynomials_cont_new(saddle_x0,saddle_x1,DF0_saddle,DF1_saddle,...
    saddle_problem0,saddle_problem1,[],Aold_saddle,Anew_saddle);
if flag <=0
    error('Could not validate saddle')
end
 
y0_point = (saddle_x0.scalar(x0.size_scalar+1:x0.size_scalar*2));
y1_point = (saddle_x1.scalar(x0.size_scalar+1:x0.size_scalar*2));
 
z0_point = (saddle_x0.scalar(x0.size_scalar*2+1:x0.size_scalar*3));
z1_point = (saddle_x1.scalar(x0.size_scalar*2+1:x0.size_scalar*3));
 
% check if validation confirmed saddle
if any(abs(y0_point)<Imin_sad) || ...
        any(abs(y1_point)<Imin_sad) || ...
        any(abs(z0_point)<Imin_sad) || ...
        any(abs(z1_point)<Imin_sad)
    bounds_too_big =1;
    %return
end
y0 = midrad(y0_point,Imin_sad);
y1 = midrad(y1_point,Imin_sad);
 
z0 = midrad(z0_point,Imin_sad);
z1 = midrad(z1_point,Imin_sad);
 
saddle_confirmed = 0;
index_saddle = find(y1.*y0<0);
if isempty(index_saddle)
    [flag_old,Imin_old]=radii_polynomials(saddle_x0,saddle_problem0,DF0_saddle,Aold_saddle);
    [flag_new,Imin_new]=radii_polynomials(saddle_x1,saddle_problem1,DF1_saddle,Anew_saddle);
    if flag_old*flag_new ==0
        error('NOOOOOOO');
    end
    I_min = min(Imin_old,Imin_new);
     
    y0 = midrad(y0_point,I_min);
    y1 = midrad(y1_point,I_min);
     
    index_saddle = find(y1.*y0<0);
    if isempty(index_saddle)
        if isempty(find(y1.*y0>0, 1))
            error('nope')
        else
            return
        end
    end
end
index_saddle= find(z0(index_saddle).*z1(index_saddle)>0);
if length(index_saddle)==1
    saddle_confirmed = 1;
    if talkative>0
        fprintf('\n       A saddle node was confirmed and is getting stored\n')
    end
elseif  isempty(index_saddle)
    saddle_confirmed = 0;
    if talkative>0
        fprintf('\n       A saddle node was not confirmed and will not be stored :(\n')
    end
else
    warning('Something is weird here: multiple variables that have a saddle node')
end
end

% function [saddle_confirmed]=validation_saddle(F_old,F_new, x0,x1,num_variable)
% global use_intlab 
% global talkative
% global refinement_saddle 
% if isempty(refinement_saddle)
%     n_intervals =300;
% else
%     n_intervals = refinement_saddle;
% end
% temp_use_intlab = use_intlab;
% use_intlab = 0;
% 
% [numerical_check] = if_saddle_numerical(F_old,F_new, x0,x1);
% 
% if ~numerical_check
%     saddle_confirmed = 0;
%     return
% end
% if ~exist('num_variable','var')|| ~exist('num_variable','var')
%     num_variable=1:x0.size_scalar;
% end
% 
% saddle_confirmed = 0;
% 
% %n_intervals= 700; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% % 
% % [numerical_check, x_prime0, x_prime1,x_prime_prime0,x_prime_prime1] = if_saddle_numerical(F_old,F_new, x0,x1);
% % step = (1/n_intervals)*(x1 - x0);
% % x_old_step = x0;
% % lin_coef_step{1} = (1/n_intervals)*(F_new.scalar_equations.linear_coef{1} - F_old.scalar_equations.linear_coef{1});
% % lin_coef_step{2} = (1/n_intervals)*(F_new.scalar_equations.linear_coef{2} - F_old.scalar_equations.linear_coef{2});
% % lin_coef_step{3} = (1/n_intervals)*(F_new.scalar_equations.linear_coef{3} - F_old.scalar_equations.linear_coef{3});
% % F_old_step = F_old;
% % numerical_check_step = zeros(n_intervals,1);
% % y0_step = zeros(n_intervals,2);
% % z0_step = zeros(n_intervals,2);
% % y1_step = zeros(n_intervals,2);
% % z1_step = zeros(n_intervals,2);
% % for i = 1: n_intervals
% %     x_new_step = x_old_step + step;
% %     F_new_step = F_old_step;
% %     
% %     F_new_step.scalar_equations.linear_coef{1} = (F_old_step.scalar_equations.linear_coef{1} + lin_coef_step{1});
% %     F_new_step.scalar_equations.linear_coef{2} = (F_old_step.scalar_equations.linear_coef{2} + lin_coef_step{2});
% %     F_new_step.scalar_equations.linear_coef{3} = (F_old_step.scalar_equations.linear_coef{3} + lin_coef_step{3});
% %     
% %     [numerical_check_step(i),x_prime0t, x_prime1t,x_prime_prime0t,x_prime_prime1t] = if_saddle_numerical(F_old_step,F_new_step, x_old_step,x_new_step);
% %     y0_step(i,:) = x_prime0t(1:2);
% %     y1_step(i,:) = x_prime1t(1:2);
% %     z0_step(i,:) = x_prime_prime0t(1:2);
% %     z1_step(i,:) = x_prime_prime1t(1:2);
% %     if numerical_check_step(i)
% %         %[~, x_prime0, x_prime1,x_prime_prime0,x_prime_prime1] = if_saddle_numerical(F_old,F_new, x0,x1);
% %         [saddle_x0, saddle_x1,y0_point,y1_point,z0_point,z1_point] = precise_computation(F_old_step, F_new_step, x_old_step,x_new_step, x_prime0, x_prime1,x_prime_prime0,x_prime_prime1);
% %         y0_step(i,:) = y0_point(1:2);
% %         y1_step(i,:) = y1_point(1:2);
% %         z0_step(i,:) = z0_point(1:2);
% %         z1_step(i,:) = z1_point(1:2);
% %     end
% %     x_old_step = x_new_step;
% %     F_old_step = F_new_step;
% % end
% % pause(1)
% 
% 
% 
% 
% 
% for i = 1: n_intervals
%     [x0_int, x1_int, F0_int, F1_int] = automatic_refinement(F_old, F_new, x0, x1, n_intervals, i); % done
%     [numerical_check, x_prime0_int, x_prime1_int,x_prime_prime0_int,...
%         x_prime_prime1_int] = if_saddle_numerical(F0_int, F1_int, x0_int, x1_int,num_variable); %done
%     
%     if ~numerical_check
%         continue
%     end
%     
%     [saddle_F0_int, saddle_F1_int] = bigger_saddle_system(F0_int, F1_int, x0_int, x1_int); %done
%     saddle_x0_int = compute_saddle(saddle_F0_int,x0_int,x_prime0_int,x_prime_prime0_int); % done
%     saddle_x1_int = compute_saddle(saddle_F1_int,x1_int,x_prime1_int,x_prime_prime1_int); % done
%     
%     % %DEBUG
%     % figure(1); hold on; plot(i,saddle_x0_int.scalar(x0.size_scalar+ num_variable),'*')
%     % figure(2); hold on; plot(i,saddle_x0_int.scalar(2*x0.size_scalar+ num_variable),'*')
%     % continue
%     % %debug
%     
%     [bool_saddle, z_crosses]=check_if_saddle_possible(saddle_x0_int,saddle_x1_int); % done
%     if bool_saddle && z_crosses
%         error('what to do now?')
%     end
%     if bool_saddle
%         use_intlab = temp_use_intlab;
%         [delta, bool_saddle_confirmed] = bound_refinement(saddle_x0_int, saddle_x1_int, saddle_F0_int, saddle_F1_int,num_variable); % done
%         use_intlab = 0;
%         if bool_saddle_confirmed
%         delta =3; bool_saddle_confirmed=0;
%         end
%         if bool_saddle_confirmed
%             use_intlab = temp_use_intlab;
%             saddle_confirmed = saddle_validation(saddle_x0_int, saddle_x1_int, saddle_F0_int, saddle_F1_int,num_variable);
%             use_intlab = temp_use_intlab;
%             if saddle_confirmed==0
%                 error('Adding intlab is problematic')
%             end
%             return
%         end
%         n_intervals_ref = ceil(1/delta)+1; %ceil(norm(norm(saddle_x0_int-saddle_x1_int))/delta)+1;
%         for j = 1: n_intervals_ref
%             [saddle_x0_small_int, saddle_x1_small_int, F0_small_int,...
%                 F1_small_int] = automatic_refinement(saddle_F0_int, ...
%                 saddle_F1_int, saddle_x0_int, saddle_x1_int, n_intervals_ref, j);% done
%             saddle_x0_small_int = Newton_2(saddle_x0_small_int,F0_small_int);
%             saddle_x1_small_int = Newton_2(saddle_x1_small_int,F1_small_int);
%             bool_saddle_small_int=check_if_saddle_possible(saddle_x0_small_int,saddle_x1_small_int); % done
%             if bool_saddle_small_int
%                 use_intlab = temp_use_intlab;
%                 bool_saddle_tot = saddle_validation(saddle_x0_small_int, saddle_x1_small_int, F0_small_int, F1_small_int,num_variable);
%                 use_intlab = 0;
%                 if bool_saddle_tot
%                     saddle_confirmed = saddle_confirmed+1;
%                     use_intlab = temp_use_intlab;
%                     return
%                 else
%                     warning('Failure to validate')
%                 end
%             end
%         end
%     else
%         if talkative >2
%             disp('False positive, not a saddle');
%         end
%     end
%     
% end
% use_intlab = temp_use_intlab;
% end
% 
% % functions to do
% % [x0_int, x1_int, F0_int, F1_int] = automatic_refinement(F_old, F_new, x0,
% % x1, n_intervals, i); % done
% % saddle_x0_int = compute_saddle(saddle_F0_int,x0_int); % done
% % [bool_saddle, z_crosses]=check_if_saddle_possible(x0_int,x1_int); % done
% % [delta, bool_saddle_confirmed] = bound_refinement(saddle_x0_int,
% % saddle_x1_int, saddle_F0_int, saddle_F1_int); % done
% % bool_saddle_tot = saddle_validation(saddle_x0_small_int, saddle_x1_small_int, F0_small_int, F1_small_int);
% 
% 
% function [x0_int, x1_int, F0_int, F1_int] = automatic_refinement(F_old, F_new, x0, x1, n_intervals, i)
% % function [x0_int, x1_int, F0_int, F1_int] = automatic_refinement(F_old, F_new, x0, x1, n_intervals, i)
% % 
% % split an interval into i subintervals and return the numerical solutions
% % and the vector fields of the subintervals
% 
% x_step = (1/n_intervals)*(x1 - x0);
% x0_int = x0 + (i-1) * x_step;
% x1_int = x0 + i * x_step;
% 
% F0_int = F_old;
% F1_int = F_old;
% 
% for j=1:3
% lin_coef_step = (1/n_intervals)*(F_new.scalar_equations.linear_coef{j} - F_old.scalar_equations.linear_coef{j});
% %lin_coef_step{2} = (1/n_intervals)*(F_new.scalar_equations.linear_coef{2} - F_old.scalar_equations.linear_coef{2});
% %lin_coef_step{3} = (1/n_intervals)*(F_new.scalar_equations.linear_coef{3} - F_old.scalar_equations.linear_coef{3});
% 
% F0_int.scalar_equations.linear_coef{j} = F_old.scalar_equations.linear_coef{j} + (i-1) *lin_coef_step;
% %F0_int.scalar_equations.linear_coef{2} = F_old.scalar_equations.linear_coef{2} + (i-1) *lin_coef_step{2};
% %F0_int.scalar_equations.linear_coef{3} = F_old.scalar_equations.linear_coef{3} + (i-1) *lin_coef_step{3};
% 
% F1_int.scalar_equations.linear_coef{j} = F_old.scalar_equations.linear_coef{j} + i *lin_coef_step;
% %F1_int.scalar_equations.linear_coef{2} = F_old.scalar_equations.linear_coef{2} + i *lin_coef_step{2};
% %F1_int.scalar_equations.linear_coef{3} = F_old.scalar_equations.linear_coef{3} + i *lin_coef_step{3};
% end
% end
% 
% function saddle_x = compute_saddle(saddle_F,x,y,z)
% % function saddle_x = compute_saddle(saddle_F,x)
% %
% % inputs: big system and small solution
% % output: big solution
% y = vec2Xi_vec(y,x);
% z = vec2Xi_vec(z,x);
% saddle_x0 = Xi_vector([x.scalar,y.scalar,z.scalar],[x.vector;y.vector;z.vector]);
% saddle_x = Newton_2(saddle_x0, saddle_F);
% end
% 
% function [bool_saddle, z_crosses]=check_if_saddle_possible(x0,x1)
% % function [bool_saddle, z_crosses]=check_if_saddle_possible(x0_int,x1_int)
% %
% % inputs: 2 big solutions
% % output: 2 bools, check if the first derivative is crossing, and if the
% % second is crossing too the zero axes
% size_scalar = x0.size_scalar /3;
% if size_scalar ~=ceil(size_scalar)
%     error('inputs off')
% end
% 
% y0 = (x0.scalar(size_scalar+1:size_scalar*2));
% y1 = (x1.scalar(size_scalar+1:size_scalar*2));
% 
% z0 = (x0.scalar(size_scalar*2+1:size_scalar*3));
% z1 = (x1.scalar(size_scalar*2+1:size_scalar*3));
% 
% index_saddle = find( y0.*y1<0);
% if ~isempty(index_saddle)
%     bool_saddle=1;
% else
%     bool_saddle=0;
%     z_crosses =0;
%     return
% end
% if any(z0(index_saddle).*z1(index_saddle)<0)
%     z_crosses=1;
% else
%     z_crosses =0;
% end
% 
% end
% 
% 
% function [delta, bool_saddle_confirmed] = bound_refinement(saddle_x0, saddle_x1, saddle_F0, saddle_F1,num_variable)
% % validate the interval and define the maximum length of the interval to be
% % able to confirm that z doesn't change sign in that interval
% bool_saddle_confirmed =0;
% 
% DF0_saddle = derivative_to_matrix(derivative(saddle_F0,saddle_x0,0));
% Aold_saddle = inv(DF0_saddle);
% DF1_saddle = derivative_to_matrix(derivative(saddle_F1,saddle_x1,0));
% Anew_saddle = inv(DF1_saddle);
% [flag,Imin,~,~,~,~,~,~,~,Ycont] = radii_polynomials_cont_new(saddle_x0,saddle_x1,DF0_saddle,DF1_saddle,...
%     saddle_F0,saddle_F1,[],Aold_saddle,Anew_saddle);
% 
% size_scalar = saddle_x0.size_scalar/3;
% z0 = (saddle_x0.scalar(size_scalar*2+1:size_scalar*3));
% z1 = (saddle_x1.scalar(size_scalar*2+1:size_scalar*3));
% z_bound = min(min(abs(z0(num_variable))),min(abs(z1(num_variable))));
% if flag>0 && Imin < z_bound
%     % validated already
%     bool_saddle_confirmed =1 ;
%     delta = 1;
%     return
%     %error('It happened!') % still to code, laziness is a bad habit
% end
% Y_extrema_max = max(Ycont(:,4));
% Y_s_max = max(Ycont(:,2));
% 
% if z_bound>Y_extrema_max%%% NEW LINE
% delta = sqrt( (z_bound - Y_extrema_max)/Y_s_max);
% else%%% NEW LINE
%     delta = 1;%%% NEW LINE
% end%%% NEW LINE
% 
% end
% 
% 
% function saddle_confirmed = saddle_validation(saddle_x0, saddle_x1, saddle_problem0, saddle_problem1,num_variable)
% 
% global talkative
% 
% DF0_saddle = derivative_to_matrix(derivative(saddle_problem0,saddle_x0,0));
% Aold_saddle = inv(DF0_saddle);
% DF1_saddle = derivative_to_matrix(derivative(saddle_problem1,saddle_x1,0));
% Anew_saddle = inv(DF1_saddle);
% [flag,Imin_sad] = radii_polynomials_cont_new(saddle_x0,saddle_x1,DF0_saddle,DF1_saddle,...
%     saddle_problem0,saddle_problem1,[],Aold_saddle,Anew_saddle);
% if flag <=0
%     error('Could not validate saddle')
% end
% size_scalar = saddle_x0.size_scalar/3;
% y0_point = (saddle_x0.scalar(size_scalar+1:size_scalar*2));
% y1_point = (saddle_x1.scalar(size_scalar+1:size_scalar*2));
% 
% z0_point = (saddle_x0.scalar(size_scalar*2+1:size_scalar*3));
% z1_point = (saddle_x1.scalar(size_scalar*2+1:size_scalar*3));
% 
% % check if validation confirmed saddle
% if any(abs(z0_point(num_variable))<Imin_sad) || ...
%         any(abs(z1_point(num_variable))<Imin_sad) ||...
%     any(abs(y0_point(num_variable))<Imin_sad) || ...
%         any(abs(y1_point(num_variable))<Imin_sad)
%     error('Bounds still too big!')
% end
% y0 = midrad(y0_point,Imin_sad);
% y1 = midrad(y1_point,Imin_sad);
% 
% z0 = midrad(z0_point,Imin_sad);
% z1 = midrad(z1_point,Imin_sad);
% 
% saddle_confirmed = 0;
% index_saddle = find(y1(num_variable).*y0(num_variable)<0);
% if isempty(index_saddle)
%     [flag_old,Imin_old]=radii_polynomials(saddle_x0,saddle_problem0,DF0_saddle,Aold_saddle);
%     [flag_new,Imin_new]=radii_polynomials(saddle_x1,saddle_problem1,DF1_saddle,Anew_saddle);
%     if flag_old*flag_new ==0
%         error('NOOOOOOO');
%     end
%     I_min = min(Imin_old,Imin_new);
%     
%     y0 = midrad(y0_point,I_min);
%     y1 = midrad(y1_point,I_min);
%     
%     index_saddle = find(y1.*y0<0);
%     if isempty(index_saddle)
%         if isempty(find(y1.*y0>0, 1))
%             error('nope')
%         else
%             return
%         end
%     end
% end
% index_saddle= find(z0(index_saddle).*z1(index_saddle)>0);
% if length(index_saddle)==1
%     saddle_confirmed = 1;
%     if talkative>0
%         fprintf('\n       A saddle node was confirmed and is getting stored\n')
%     end
% elseif  isempty(index_saddle)
%     saddle_confirmed = 0;
%     if talkative>0
%         fprintf('\n       A saddle node was not confirmed and will not be stored :(\n')
%     end
% else
%     warning('Something is weird here: multiple variables that have a saddle node')
% end
% end
% 
% %%%%%%%%%%%% NOT USED ANYMORE
% function [saddle_confirmed]=validation_saddle_nope(F_old,F_new, x0,x1)
% % function [saddle_confirmed]=validation_saddle(F_old,F_new, x0,x1)
% %
% % validation of a saddle node between x0 and x1
% % outline:
% % numerical check
% % construction of bigger system
% % Newton
% % validation of bigger system
% % end checks
% % storage
% %
% % INPUT
% % F_old, F_new      full_problems, square problems, including continuation
% %                   equation
% % x0, x1            Xi_vectors, validated solutions to F_old and F_new
% %                   respectively
% % OUTPUT
% % saddle_confirmed  bool, 1 if the saddle bifurcation was validated, 0
% %                   otherwise
% 
% % other outputs to add
% % saddle_index_var saddle_x0 saddle_x1 y0 y1 z0 z1
% % %
% % % % new pseudocode:
% % % % numerical check
% % % % test on [x0,x1]
% % % % if radmin > abs(y) or abs(z) AND || x0  - x1 || > delta
% % % %   if numerical check (x0, midpoint)
% % % %    test on [x0, mid_point]
% % % %   else if numerical check (midpoint, x1)
% % % %    test on [mid_point, x0]
% % % %   else
% % % %     wtf
% % % %    end
% % % % else
% % % %    impossible validation
% % % % end
% 
% global talkative
% 
% % numerical check
% [numerical_check, x_prime0, x_prime1,x_prime_prime0,x_prime_prime1] = if_saddle_numerical(F_old,F_new, x0,x1);
% 
% if ~numerical_check
%     saddle_confirmed = 0;
%     return
% end
% % [saddle_confirmed,bounds_too_big] = actual_validation(F_old, F_new, x0,x1, x_prime0, x_prime1,x_prime_prime0,x_prime_prime1);
% 
% 
% % % if saddle_confirmed ==1
% % %     return
% % % elseif bounds_too_big ==1 && norm(norm(x0 - x1))>10^-9
% % %     % divide segment in half
% %     x_mid = 0.5*(x0+x1);
% %     F_mid = F_old;
% %     F_mid.scalar_equations.linear_coef{1} = 0.5*(F_old.scalar_equations.linear_coef{1} + F_new.scalar_equations.linear_coef{1});
% %     F_mid.scalar_equations.linear_coef{2} = 0.5*(F_old.scalar_equations.linear_coef{2} + F_new.scalar_equations.linear_coef{2});
% %     F_mid.scalar_equations.linear_coef{3} = 0.5*(F_old.scalar_equations.linear_coef{3} + F_new.scalar_equations.linear_coef{3});
% %
% %     [numerical_check1] = if_saddle_numerical(F_old,F_mid, x0,x_mid);
% %     [numerical_check2] = if_saddle_numerical(F_mid,F_new, x_mid,x1);
% %     if numerical_check1 ==0 && numerical_check2 ==0
% %         error('something went weird');
% %     elseif numerical_check1
% %         [saddle_confirmed] = actual_validation(F_old, F_mid, x0,x_mid, x_prime0, x_prime1,x_prime_prime0,x_prime_prime1);
% %         %[saddle_confirmed]=validation_saddle(F_old,F_mid, x0,x_mid);
% %     else
% %         [saddle_confirmed] = actual_validation(F_mid, F_new, x_mid,x1, x_prime0, x_prime1,x_prime_prime0,x_prime_prime1);
% %         %[saddle_confirmed]=validation_saddle(F_mid,F_new, x_mid,x1);
% %     end
% % % end
% z0 = x_prime_prime0(1:x0.size_scalar);
% z1 = x_prime_prime1(1:x0.size_scalar);
% bound = min( min ( abs( z0), abs(z1)))/2;
% 
% %maximum_length = test_validation(bound,F_old, F_new, x0,x1, x_prime0, x_prime1,x_prime_prime0,x_prime_prime1);
% 
% n_intervals = 30;%ceil(norm(norm(x0-x1))/maximum_length)+1;
% 
% step = (1/n_intervals)*(x1 - x0);
% x_old_step = x0;
% lin_coef_step{1} = (1/n_intervals)*(F_new.scalar_equations.linear_coef{1} - F_old.scalar_equations.linear_coef{1});
% lin_coef_step{2} = (1/n_intervals)*(F_new.scalar_equations.linear_coef{2} - F_old.scalar_equations.linear_coef{2});
% lin_coef_step{3} = (1/n_intervals)*(F_new.scalar_equations.linear_coef{3} - F_old.scalar_equations.linear_coef{3});
% F_old_step = F_old;
% numerical_check_step = zeros(n_intervals,1);
% y0_step = zeros(n_intervals,2);
% z0_step = zeros(n_intervals,2);
% y1_step = zeros(n_intervals,2);
% z1_step = zeros(n_intervals,2);
% for i = 1: n_intervals
%     x_new_step = x_old_step + step;
%     F_new_step = F_old_step;
%     
%     F_new_step.scalar_equations.linear_coef{1} = (F_old_step.scalar_equations.linear_coef{1} + lin_coef_step{1});
%     F_new_step.scalar_equations.linear_coef{2} = (F_old_step.scalar_equations.linear_coef{2} + lin_coef_step{2});
%     F_new_step.scalar_equations.linear_coef{3} = (F_old_step.scalar_equations.linear_coef{3} + lin_coef_step{3});
%     
%     [numerical_check_step(i),x_prime0t, x_prime1t,x_prime_prime0t,x_prime_prime1t] = if_saddle_numerical(F_old_step,F_new_step, x_old_step,x_new_step);
%     y0_step(i,:) = x_prime0t(1:2);
%     y1_step(i,:) = x_prime1t(1:2);
%     z0_step(i,:) = x_prime_prime0t(1:2);
%     z1_step(i,:) = x_prime_prime1t(1:2);
%     if numerical_check_step(i)
%         %[~, x_prime0, x_prime1,x_prime_prime0,x_prime_prime1] = if_saddle_numerical(F_old,F_new, x0,x1);
%         [saddle_x0, saddle_x1,y0_point,y1_point,z0_point,z1_point] = precise_computation(F_old_step, F_new_step, x_old_step,x_new_step, x_prime0, x_prime1,x_prime_prime0,x_prime_prime1);
%         y0_step(i,:) = y0_point(1:2);
%         y1_step(i,:) = y1_point(1:2);
%         z0_step(i,:) = z0_point(1:2);
%         z1_step(i,:) = z1_point(1:2);
%     end
%     x_old_step = x_new_step;
%     F_old_step = F_new_step;
% end
% pause(1)
% end
% 
% function maximum_step = test_validation(bound,F_old, F_new, x0,x1, x_prime0, x_prime1,x_prime_prime0,x_prime_prime1)
% 
% [saddle_problem0,saddle_problem1] = bigger_saddle_system(F_old,F_new,x0, x1);
% 
% y0 = vec2Xi_vec(x_prime0,x0);
% y1 = vec2Xi_vec(x_prime1,x1);
% z0 = vec2Xi_vec(x_prime_prime0,x0);
% z1 = vec2Xi_vec(x_prime_prime1,x1);
% saddle_x0 = Xi_vector([x0.scalar,y0.scalar,z0.scalar],[x0.vector;y0.vector;z0.vector]);
% saddle_x1 = Xi_vector([x1.scalar,y1.scalar,z1.scalar],[x1.vector;y1.vector;z1.vector]);
% 
% saddle_x0=Newton_2(saddle_x0,saddle_problem0);
% saddle_x1=Newton_2(saddle_x1,saddle_problem1);
% 
% % reset x to initial values, it cannot be changed from what it was
% saddle_x0.scalar(1:x0.size_scalar) = x0.scalar;
% saddle_x0.vector(1:x0.size_vector,:) = x0.vector;
% saddle_x1.scalar(1:x1.size_scalar) = x1.scalar;
% saddle_x1.vector(1:x1.size_vector,:) = x1.vector;
% 
% % VALIDATION
% 
% DF0_saddle = derivative_to_matrix(derivative(saddle_problem0,saddle_x0,0));
% Aold_saddle = inv(DF0_saddle);
% DF1_saddle = derivative_to_matrix(derivative(saddle_problem1,saddle_x1,0));
% Anew_saddle = inv(DF1_saddle);
% [flag,Imin,~,~,~,~,~,~,~,Ycont] = radii_polynomials_cont_new(saddle_x0,saddle_x1,DF0_saddle,DF1_saddle,...
%     saddle_problem0,saddle_problem1,[],Aold_saddle,Anew_saddle);
% if flag>0 && Imin < bound
%     % validated already
%     pause(1)
% end
% Y_extrema_max = max(Ycont(:,4));
% Y_s_max = max(Ycont(:,2));
% 
% maximum_step = sqrt( (bound - Y_extrema_max)/Y_s_max);
% 
% end
% 
% function [saddle_x0, saddle_x1,y0_point,y1_point,z0_point,z1_point] = precise_computation(F_old, F_new, x0,x1, x_prime0, x_prime1,x_prime_prime0,x_prime_prime1)
% % outline
% % construction of bigger system
% % Newton
% % validation of bigger system
% % end checks
% 
% saddle_confirmed = 0;
% bounds_too_big = 0;
% 
% % setting for validation
% % construct bigger system
% [saddle_problem0,saddle_problem1, DxxHyy] = bigger_saddle_system(F_old,F_new,x0, x1); % DxxHyy still to debug: dimensions of linear coefficients not matching
% 
% % not in use anymore TEMPORARY
% % index_old_x = 2; % index_old_x is the index of the scalar variable w.r.t. which we have the bifurcation
% 
% y0 = vec2Xi_vec(x_prime0,x0);
% y1 = vec2Xi_vec(x_prime1,x1);
% z0 = vec2Xi_vec(x_prime_prime0,x0);
% z1 = vec2Xi_vec(x_prime_prime1,x1);
% saddle_x0 = Xi_vector([x0.scalar,y0.scalar,z0.scalar],[x0.vector;y0.vector;z0.vector]);
% saddle_x1 = Xi_vector([x1.scalar,y1.scalar,z1.scalar],[x1.vector;y1.vector;z1.vector]);
% 
% % solve for z variables
% % by direct computation - not necessary anymore
% % z0 = second_archlength_der(F_old,F_old,F_new,DxxHyy,x0,y0,saddle_x0);
% % z1 = second_archlength_der(F_new,F_old,F_new,DxxHyy,x1,y1,saddle_x1);
% %
% % z0 = vec2Xi_vec(z0,x0);
% % z1 = vec2Xi_vec(z1,x1);
% %
% % saddle_x0 = Xi_vector([x0.scalar,y0.scalar,z0.scalar],[x0.vector;y0.vector;z0.vector]);
% % saddle_x1 = Xi_vector([x1.scalar,y1.scalar,z1.scalar],[x1.vector;y1.vector;z1.vector]);
% 
% saddle_x0=Newton_2(saddle_x0,saddle_problem0);
% saddle_x1=Newton_2(saddle_x1,saddle_problem1);
% y0_point = (saddle_x0.scalar(x0.size_scalar+1:x0.size_scalar*2));
% y1_point = (saddle_x1.scalar(x0.size_scalar+1:x0.size_scalar*2));
% 
% z0_point = (saddle_x0.scalar(x0.size_scalar*2+1:x0.size_scalar*3));
% z1_point = (saddle_x1.scalar(x0.size_scalar*2+1:x0.size_scalar*3));
% end
% 
% 
% function [saddle_confirmed, bounds_too_big] = actual_validation(F_old, F_new, x0,x1, x_prime0, x_prime1,x_prime_prime0,x_prime_prime1)
% % outline
% % construction of bigger system
% % Newton
% % validation of bigger system
% % end checks
% 
% saddle_confirmed = 0;
% bounds_too_big = 0;
% 
% % setting for validation
% % construct bigger system
% [saddle_problem0,saddle_problem1, DxxHyy] = bigger_saddle_system(F_old,F_new,x0, x1); % DxxHyy still to debug: dimensions of linear coefficients not matching
% 
% % not in use anymore TEMPORARY
% % index_old_x = 2; % index_old_x is the index of the scalar variable w.r.t. which we have the bifurcation
% 
% y0 = vec2Xi_vec(x_prime0,x0);
% y1 = vec2Xi_vec(x_prime1,x1);
% z0 = vec2Xi_vec(x_prime_prime0,x0);
% z1 = vec2Xi_vec(x_prime_prime1,x1);
% saddle_x0 = Xi_vector([x0.scalar,y0.scalar,z0.scalar],[x0.vector;y0.vector;z0.vector]);
% saddle_x1 = Xi_vector([x1.scalar,y1.scalar,z1.scalar],[x1.vector;y1.vector;z1.vector]);
% 
% % solve for z variables
% % by direct computation - not necessary anymore
% % z0 = second_archlength_der(F_old,F_old,F_new,DxxHyy,x0,y0,saddle_x0);
% % z1 = second_archlength_der(F_new,F_old,F_new,DxxHyy,x1,y1,saddle_x1);
% %
% % z0 = vec2Xi_vec(z0,x0);
% % z1 = vec2Xi_vec(z1,x1);
% %
% % saddle_x0 = Xi_vector([x0.scalar,y0.scalar,z0.scalar],[x0.vector;y0.vector;z0.vector]);
% % saddle_x1 = Xi_vector([x1.scalar,y1.scalar,z1.scalar],[x1.vector;y1.vector;z1.vector]);
% 
% saddle_x0=Newton_2(saddle_x0,saddle_problem0);
% saddle_x1=Newton_2(saddle_x1,saddle_problem1);
% 
% % reset x to initial values, it cannot be changed from what it was
% saddle_x0.scalar(1:x0.size_scalar) = x0.scalar;
% saddle_x0.vector(1:x0.size_vector,:) = x0.vector;
% saddle_x1.scalar(1:x1.size_scalar) = x1.scalar;
% saddle_x1.vector(1:x1.size_vector,:) = x1.vector;
% 
% % VALIDATION
% 
% % % with a particular Newton method
% % iter_Newton = 0;
% % while iter_Newton < 30 && max(norm(apply(saddle_problem0,saddle_x0)))>min_res_N
% %     saddle_x0 = saddle_x0 - vec2Xi_vec( derivative_to_matrix(derivative(saddle_problem0,saddle_x0,0))\ Xi_vec2vec(apply(saddle_problem0,saddle_x0)) , saddle_x0);
% %     % saddle_x0_vec = Xi_vec2vec(saddle_x0);
% %     % saddle_x = saddle_x0_vec - derivative_to_matrix(derivative(saddle_problem0,saddle_x0,0))\ Xi_vec2vec(apply(saddle_problem0,saddle_x0));
% %     % saddle_x(index_old_x) = saddle_x0_vec(index_old_x);
% %     % saddle_x0 = vec2Xi_vec(saddle_x,saddle_x0);
% % end
% % if max(norm(apply(saddle_problem0,saddle_x0)))>min_res_N
% %     error('Newton for the saddle node validation at x0 did not converge')
% % end
% % iter_Newton = 0;
% % while iter_Newton < 30 && max(norm(apply(saddle_problem1,saddle_x1)))>min_res_N
% %     saddle_x1_vec = Xi_vec2vec(saddle_x1);
% %     saddle_x = saddle_x1_vec - derivative_to_matrix(derivative(saddle_problem1,saddle_x1,0))\ Xi_vec2vec(apply(saddle_problem1,saddle_x1));
% %     saddle_x(index_old_x) = saddle_x1_vec(index_old_x);
% %     saddle_x1 = vec2Xi_vec(saddle_x,saddle_x1);
% % end
% % if max(norm(apply(saddle_problem1,saddle_x1)))>min_res_N
% %     error('Newton for the saddle node validation at x1 did not converge')
% % end
% 
% global talkative
% 
% % validate segment
% 
% DF0_saddle = derivative_to_matrix(derivative(saddle_problem0,saddle_x0,0));
% Aold_saddle = inv(DF0_saddle);
% DF1_saddle = derivative_to_matrix(derivative(saddle_problem1,saddle_x1,0));
% Anew_saddle = inv(DF1_saddle);
% [flag,Imin_sad] = radii_polynomials_cont_new(saddle_x0,saddle_x1,DF0_saddle,DF1_saddle,...
%     saddle_problem0,saddle_problem1,[],Aold_saddle,Anew_saddle);
% if flag <=0
%     error('Could not validate saddle')
% end
% 
% y0_point = (saddle_x0.scalar(x0.size_scalar+1:x0.size_scalar*2));
% y1_point = (saddle_x1.scalar(x0.size_scalar+1:x0.size_scalar*2));
% 
% z0_point = (saddle_x0.scalar(x0.size_scalar*2+1:x0.size_scalar*3));
% z1_point = (saddle_x1.scalar(x0.size_scalar*2+1:x0.size_scalar*3));
% 
% % check if validation confirmed saddle
% if any(abs(y0_point)<Imin_sad) || ...
%         any(abs(y1_point)<Imin_sad) || ...
%         any(abs(z0_point)<Imin_sad) || ...
%         any(abs(z1_point)<Imin_sad)
%     bounds_too_big =1;
%     %return
% end
% y0 = midrad(y0_point,Imin_sad);
% y1 = midrad(y1_point,Imin_sad);
% 
% z0 = midrad(z0_point,Imin_sad);
% z1 = midrad(z1_point,Imin_sad);
% 
% saddle_confirmed = 0;
% index_saddle = find(y1.*y0<0);
% if isempty(index_saddle)
%     [flag_old,Imin_old]=radii_polynomials(saddle_x0,saddle_problem0,DF0_saddle,Aold_saddle);
%     [flag_new,Imin_new]=radii_polynomials(saddle_x1,saddle_problem1,DF1_saddle,Anew_saddle);
%     if flag_old*flag_new ==0
%         error('NOOOOOOO');
%     end
%     I_min = min(Imin_old,Imin_new);
%     
%     y0 = midrad(y0_point,I_min);
%     y1 = midrad(y1_point,I_min);
%     
%     index_saddle = find(y1.*y0<0);
%     if isempty(index_saddle)
%         if isempty(find(y1.*y0>0, 1))
%             error('nope')
%         else
%             return
%         end
%     end
% end
% index_saddle= find(z0(index_saddle).*z1(index_saddle)>0);
% if length(index_saddle)==1
%     saddle_confirmed = 1;
%     if talkative>0
%         fprintf('\n       A saddle node was confirmed and is getting stored\n')
%     end
% elseif  isempty(index_saddle)
%     saddle_confirmed = 0;
%     if talkative>0
%         fprintf('\n       A saddle node was not confirmed and will not be stored :(\n')
%     end
% else
%     warning('Something is weird here: multiple variables that have a saddle node')
% end
% end