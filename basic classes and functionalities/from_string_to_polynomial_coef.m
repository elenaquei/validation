function vector_field = from_string_to_polynomial_coef(string_s, number_scalars, number_function_unknowns)
% function vector_field = from_string_to_polynomial_coef(string_s, number_scalars, number_functions_unknowns)
%
% INPUT
% string_s             a string with the vector field
% number_scalars       integer, number of scalar unknows (DEFAULT: the
%                      highest li)    
% number_functions     integer, number of time-dependent unknowns (DEFAULT: the
%                      highest xi)    
% OUTPUT
% vector_filed         instance of the class polynomial_coef
% disp(vector_field)   for checking
%
% string_s must be something of the form
%  + 3 dot x2 + x2 x3 l1 l2
%  + dot x1 + l1 l2 x2 x3
%  + 45 x1 - 2 dot x1
% 
% the scalar unknowns are l1, l2, ---- and the time-dependent solutions are
% x1, x2, --- and must be used consequently if the two optional parameters
% are not inserted
% every line needs to start with a sign (+ or -), 1s are not necessary, the
% derivative with respect to time is signaled as 'dot xi', there must be
% spaces inbetween every element. Underscores are not necessary, product
% signs are not necessary, i.e. the following two lines are equivalent:
% + 3 dot x2 + x2^4 x3 l1 l2
% and
% 3 dot x_2 + x2 ^ 4 x3 * l_1 * l_2  (mixed notation possible)
% NO PARENTHESIS ALLOWED

% struncture: 
% 0) delete * _ 
% 1) check how many lines there are
% 2) if number_scalars and/or number_function_unknowns are not inserted,
% search for them (M and N)
% 0b) check that no other element is inserted than
% 1,2,---0,x,l,d,o,t,+,-,space,\,n,^
% 3) look at every term and decript it
%      3a) get value{i}(j)
%      3b) get any dot, store it and delete it
%      3c) get l1....lm, store it and delete it
%      3d) get x1....xn, store it and delete it
% 4) put all together
% 5) double-check: print the vector field
global talkative

% delete unwanted signs 
index_star = find(string_s =='*'); 
string_s(index_star) = [];
index_star = find(string_s =='_');
string_s(index_star) = [];

index_space = find(string_s ==' '); % delete all spacing 
string_s(index_space) = [];

% if there is any symbol that is not 1,2,---0,x,l,d,o,t,+,-,space,\,n, then crash
% here
if any(((string_s == '1') + (string_s=='2') + (string_s=='3') + (string_s=='4') ...
        + (string_s=='5') + (string_s=='6') + (string_s=='7') + (string_s=='8')...
         + (string_s=='9') + (string_s=='0') + (string_s=='x') + (string_s=='l')...
          + (string_s=='d') + (string_s=='o') + (string_s=='t') + (string_s=='+')...
           + (string_s=='-') + (string_s=='\') + (string_s=='n')+ (string_s=='^')+ (string_s=='.'))~=1)
    error('Unrecognised character')
end

if length(find(string_s=='d'))~=length(strfind(string_s,'dotx'))
    error('Unrecognised sequence of characters')
end

%  number_scalars, number_function_unknowns setted 
if nargin<2 || isempty(number_scalars)
    index_l = find(string_s=='l');
    number_scalars = 1;
    for i =1:length(index_l)
        if isempty(str2num(string_s(index_l(i)+1)))
            error('No scalar numbering')
        end
        number_scalars = max(number_scalars,str2num(string_s(index_l(i)+1)));
    end
end
if nargin<3 || isempty(number_function_unknowns)
    index_x = find(string_s=='x');
    number_function_unknowns = 1;
    for i =1:length(index_x)
        if isempty(str2num(string_s(index_x(i)+1)))
            error('Index of vector unknown not defined in %d', index_x(i))
        end
        [~,number_of_x] = length_double(string_s,index_x(i)+1);
        number_function_unknowns = max(number_function_unknowns,number_of_x);%str2num(string_s(index_x(i)+1)));
    end
end

% find out how many equations there are
index_new_line = find((string_s(1:end-1) =='/') + (string_s(2:end)=='n'));
index_new_line = [1,index_new_line];
if string_s(end) ~= 'n'
    index_new_line = [index_new_line, length(string_s)];
end
num_equations = length(index_new_line)-1;
string_eq = cell(num_equations,1);
for i = 1:num_equations
    string_eq{i} = string_s(index_new_line(i) : index_new_line(i+1));
    index_unused = find((string_eq{i}=='\')+(string_eq{i}=='n'));
    string_eq{i}(index_unused) =[];
end

% find out how many terms there are 
n_terms = zeros(1,num_equations);
start_terms = cell(1,num_equations);
end_terms = cell(1,num_equations);


% preallocation of elementsvalue = cell(num_equations,1);
dot = cell(num_equations,1);
power_vector = cell(num_equations,1);
delay = cell(num_equations,1);
power_scalar = cell(num_equations,1);
value= cell(num_equations,1);
    % value % cell{n_equations}(n_terms) -- coefficient of the term
    % dot % cell{n_equations}{n_terms}(variables,:) -- number of derivatives applied to the given variable 
    % power_vector % cell{n_equations}{n_terms}(variables,:)
    % delay % cell{n_equations}{n_terms}(variables,:)
    % power_scalar % cell{n_equations}(scalar_variables,n_terms)


for i = 1:num_equations
    if (string_eq{i}(1)~='+') && (string_eq{i}(1) ~='-')
        string_eq{i}= sprintf('+%s',string_eq{i});
    end
    position_signs = find((string_eq{i}=='+')+ (string_eq{i}=='-'));
    n_terms(i) = length(position_signs);
    % where they are
    start_terms{i} = position_signs;
    end_terms{i} = [position_signs(2:end)-1, length(string_eq{i})];
    
    value{i} = zeros(n_terms(i),1);
    dot{i} =  cell(n_terms(i),1);%(variables,1)
    delay{i} =  cell(n_terms(i),1);
    power_vector{i} =  cell(n_terms(i),1);% (variables,1)
    power_scalar{i}=zeros(number_scalars,n_terms(i));
    
    for j = 1:n_terms(i)
        % preallocation of elements
        dot{i}{j} =  zeros(number_function_unknowns,1);
        delay{i}{j} =  zeros(number_function_unknowns,1);
        power_vector{i}{j} =  zeros (number_function_unknowns,1);
        
        % and what they are
        string_term = string_eq{i}(start_terms{i}(j):end_terms{i}(j));
        % determine the scalar coefficient in front
        if isempty(str2num(string_term(1:2)))
            value{i}(j) = str2num(strcat(string_term(1),'1'));
            string_term = string_term(2:end);
        else
            % index_value = 2;
            % while  index_value+1<=length(string_term) && ~isempty(str2num(string_term(1:index_value+1)))
            %     index_value = index_value+1;
            % end
            [index_end,double] = length_double(string_term, 1);
            value{i}(j) = double;%str2num(string_term(1:index_value));
            string_term = string_term(index_end+1:end);
        end
        % go along and find all elements
        % each time somethign is found, it is also deleted
        while ~isempty(string_term)
           % dots
           %     each time with following powers
           index_dot = strfind(string_term,'dotx');
           if ~isempty(index_dot)
               if length(index_dot)>1
                   warning('At the moment, the code is still not ready for multiple dots, but this function will compile')
               end
               index_dot = index_dot(1);
               [end_number_x,number_x] = length_double(string_term, index_dot+4);%str2num(string_term(index_dot+4));
               if isempty(number_x)
                   error('Dot of unknown element at %s',string_term(index_dot+(0:3)));
               end
               
               % length_dot = 4;
               %     each time with following powers
               if string_s(end_number_x+1)=='^'
                   power_starts = end_number_x+2;
                   if isempty(str2num(string_term(power_starts)))
                       error('Signal of ^ but no digit afterwards');
                   end
                   %length_dot = 6;
                   %while  power_starts+1<=length(string_term) && ~isempty(str2num(string_term(index_dot+6:power_starts+1)))
                   %    power_starts = power_starts+1;
                   %    length_dot = length_dot+1;
                   %end
                   [end_dot,dot{i}{j}(number_x)] = length_double(string_term,power_starts);%str2num(string_term(index_dot+6:power_starts));% % MISSING : change power storage
               else
                   dot{i}{j}(number_x) =1;% % MISSING : change power storage
                   end_dot = end_number_x;
               end
               % each time deleating the used parts
               string_term (index_dot:end_dot) =[];
               continue % we want to  get rid of all the dots before going on with xs (otherwise trouble)
           end
           
           % xs
           %     each time with following powers
           index_x = strfind(string_term,'x');
           if ~isempty(index_x)
               if length(index_x)>1
                   index_x = index_x(1);
               end
               [end_x,number_x] = length_double(string_term, index_x+1);%str2num(string_term(index_x+1));
               if isempty(number_x)
                   error('Dot of unknown element at %s',string_term(index_x+(0:1)));
               end
               
               % length_x = 1;
               %     each time with following powers
               if length(string_term)>=end_x+1 && string_term(end_x+1)=='^'
                   power_starts = end_x+2;
                   if isempty(str2num(string_term(power_starts)))
                       error('Signal of ^ but no digit afterwards');
                   end
                   %length_x = 3;
                   %while  power_starts+1<=length(string_term) && ~isempty(str2num(string_term(index_x+3:power_starts+1)))
                   %    power_starts = power_starts+1;
                   %    length_x = length_x+1;
                   %end
                   [end_x,power_vector{i}{j}(number_x)] = length_double(string_term,power_starts);
                   %str2num(string_term(index_x+3:power_starts));% % MISSING : change power storage
               else
                   power_vector{i}{j}(number_x) =1;% % MISSING : change power storage
                   %end_x = index_x+1;
               end
               % each time deleating the used parts
               string_term (index_x:end_x) =[];
           end
           
           % ls
           %     each time with following powers
%            index_l = strfind(string_term,'l');
%            if ~isempty(index_l)
%                if length(index_l)>1
%                    index_l = index_l(1);
%                end
%                number_l = str2num(string_term(index_l+1));
%                if isempty(number_l)
%                    error('Dot of unknown element at %s',string_term(index_l+(0:1)));
%                end
%                length_l = 1;
%                %     each time with following powers
%                if length(string_term)>=index_l+2 &&string_term(index_l+2)=='^'
%                    power_starts = index_l+3;
%                    if isempty(str2num(string_term(index_l+3)))
%                        error('Signal of ^ but no digit afterwards');
%                    end
%                    %length_l = 3;
%                    %while  power_starts+1<=length(string_term) && ~isempty(str2num(string_term(index_l+3:power_starts+1)))
%                    %    power_starts = power_starts+1;
%                    %    length_l = length_l+1;
%                    %end
%                    [end_l,power_scalar{i}(number_l,j)] = length_double(string_term,power_starts);
%                    %str2num(string_term(index_l+3:power_starts));
%                else
%                    power_scalar{i}(number_l,j) = 1;
%                    end_l = index_l + 1;
%                end
%                % each time deleating the used parts
%                string_term (index_l:end_l) =[];
%            end
           index_l = strfind(string_term,'l');
           if ~isempty(index_l)
               if length(index_l)>1
                   index_l = index_l(1);
               end
               [end_l,number_l] = length_double(string_term, index_l+1);%str2num(string_term(index_x+1));
               if isempty(number_l)
                   error('Dot of unknown element at %s',string_term(index_l+(0:1)));
               end
               
               % length_x = 1;
               %     each time with following powers
               if length(string_term)>=end_l+1 && string_term(end_l+1)=='^'
                   power_starts = end_l+2;
                   if isempty(str2num(string_term(power_starts)))
                       error('Signal of ^ but no digit afterwards');
                   end
                   %length_x = 3;
                   %while  power_starts+1<=length(string_term) && ~isempty(str2num(string_term(index_x+3:power_starts+1)))
                   %    power_starts = power_starts+1;
                   %    length_x = length_x+1;
                   %end
                   [end_l,power_scalar{i}(number_l,j)] = length_double(string_term,power_starts);
                   %str2num(string_term(index_x+3:power_starts));% % MISSING : change power storage
               else
                   power_scalar{i}(number_l,j) = 1;
                   %end_x = index_x+1;
               end
               % each time deleating the used parts
               string_term (index_l:end_l) =[];
           end
           
        end
    end
end

vector_field = polynomial_coefs(number_scalars, number_function_unknowns, num_equations, ...
                n_terms, value,power_scalar,power_vector, dot, delay);
if talkative>0
    disp(vector_field);
end
end

function [index_end,double] = length_double(string, index)
% function [index_end,double] = length_double(string, index)
%
% this function checks how long is the numeral starting at index in string
% and returns the index of the last element of the numeral

length_string = length(string);
index_start = index;
sign = 1;

if strcmp(string(index_start),'+')
    index_start = index_start+1;
elseif strcmp(string(index),'-')
    sign = -1;
    index_start = index_start+1;
end

if strcmp(string(index_start),'.')
    index_start = index_start+1;
end

if isempty(str2num(string(index:length_string))) || isnan(str2double(string(index:length_string)))
    for i = index_start:length_string
        if isempty(str2num(string(index:i)))
            break
        end
    end
    index_end = i-1;
else
    index_end = length_string;
end

double = str2double (string(index:index_end));

end