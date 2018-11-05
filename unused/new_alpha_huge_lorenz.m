function [alpha_cont,s]=new_alpha_huge_lorenz(DIM,s,nodes)
% function [alpha_cont,s]=new_alpha_huge_lorenz(DIM,s,nodes)
%
% INPUT
% DIM      number of correlated Lorenz systems
% s        string of name of saved file
% nodes    number of nodes
% OUTPUT
% s        name of the saved system
%
% alpha for super_huge Lorenz:
% DIM correlated Lorenz systems
% 
% (x1,y1,z1)...(xDIM,yDIM,zDIM)
%
% with intertwined dependency

% working as of 9th July 2019

if nargin<3 || isempty(s)
    s=sprintf('huge_lorenz%d',DIM);
end

num_vec=DIM*3;
sigma=10;beta=8/3; %pho=28;

string_lorenz = '- dot x_3i+1 + sigma l1 x_3i+2 - sigma l1 x_3i+1 \n - dot x_3i+2 + pho l1 x_3i+1 - l1 x_3i+1 x_3i+3 - l1 x_3i+5 \n - dot x_3i+3 + l1 x_3i+1 x_3i+2 - beta l1 x_3i+3\n'; % general lorenz
string_lorenz_vars = strrep(string_lorenz, 'sigma' , num2str(sigma)); % plugging in sigma
string_lorenz_vars = strrep(string_lorenz_vars, 'beta' , num2str(beta)); % plugging in beta
%string_lorenz_pho = strrep(string_lorenz_vars, 'pho', num2str(pho_null)); % for point wise first system, plugging in pho
string_lorenz_cont = strrep(string_lorenz_vars, 'pho', 'l2'); % setting pho as the second scalar variable

string_vf='';
for i=1:DIM
    string_lorenz_i=strrep(string_lorenz_cont, '_3i+1', num2str(mod(3*(i-1)+1-1,3*DIM)+1));
    string_lorenz_i=strrep(string_lorenz_i, '_3i+2', num2str(mod(3*(i-1)+2-1,3*DIM)+1));
    string_lorenz_i=strrep(string_lorenz_i, '_3i+3', num2str(mod(3*(i-1)+3-1,3*DIM)+1));
    string_lorenz_i=strrep(string_lorenz_i, '_3i+5', num2str(mod(3*(i-1)+5-1,3*DIM)+1));
    string_vf=strcat(string_vf,string_lorenz_i);
end

alpha_cont = from_string_to_polynomial_coef(string_vf);

return

%%
if nargin==1
    epsilon=10^-4;
end

n_term = repmat([3,5,3],1,DIM);
power_scal = cell(num_vec,1);
power_vector = cell(num_vec,1);
value = cell(num_vec,1);%{n_equations}(n_terms)
dot = cell(num_vec,1);
delay = dot;
for i = 1:num_vec
    dot{i} = cell(n_term(i),1);
    delay{i} = dot{i};
    for j = 1:n_term(i)
        dot{i}{j} = zeros(num_vec,1);
    end
    dot{i}{1}(i) = 1;
    power_vector{i} = cell(n_term(i),1);
    power_scal{i} = zeros(2,n_term(i));
    power_scal{i}(1,2:end)=1;
end
x = [1,0,0].';
y = [0,1,0].';
z = [0,0,1].';
for i = 0:DIM-1
    power_scal{i*3+2}(2,5)=1;
    
    power_vector{i*3+1}{1} = zeros(num_vec,1);
    power_vector{i*3+1}{2} = zeros(num_vec,1);
    power_vector{i*3+1}{3} = zeros(num_vec,1);
    
    power_vector{i*3+2}{1} = zeros(num_vec,1);
    power_vector{i*3+2}{2} = zeros(num_vec,1);
    power_vector{i*3+2}{3} = zeros(num_vec,1);
    power_vector{i*3+2}{4} = zeros(num_vec,1);
    power_vector{i*3+2}{5} = zeros(num_vec,1);
    
    power_vector{i*3+3}{1} = zeros(num_vec,1);
    power_vector{i*3+3}{2} = zeros(num_vec,1);
    power_vector{i*3+3}{3} = zeros(num_vec,1);
    
    coord = i*3+(1:3);
    power_vector{i*3+1}{2}(coord) = y;
    power_vector{i*3+1}{3}(coord) = x;
    
    power_vector{i*3+2}{2}(coord) = x;
    power_vector{i*3+2}{3}(coord) = x+z;
    power_vector{i*3+2}{4}(coord) = y;
    coord_4 =mod(coord+3,num_vec+1);
    if coord_4(1) ==0
        coord_4 = coord_4+1;
    end
    power_vector{i*3+2}{5}(coord_4) = 2*y;
    
    power_vector{i*3+3}{2}(coord) = x+y;
    power_vector{i*3+3}{3}(coord) = z;
    
    value{i*3+1} = [1,-sigma,sigma];
    value{i*3+2} = [1,-pho,1,1,1];
    value{i*3+3} = [1,-1,beta];
end

alpha_vec = polynomial_coefs(2, num_vec, num_vec, ...
    n_term, value,power_scal,power_vector, dot);

pol = polynomial_coefs(2,num_vec,0,zeros(0),cell(0),cell(0),cell(0),cell(0));
lin_coef = cell(3,1);
lin_coef{1} = ones(1,2);
lin_coef{2} = ones(1,num_vec,nodes*2+1);
lin_coef{3} = [0];
alpha_scal = scalar_eq(1,0,2,num_vec,lin_coef,pol);

alpha_cont = full_problem(alpha_scal, alpha_vec);

if nargout ==2
    save(s,'alpha_cont');
end


% value=cell(num_vec,1);
% powers_scalar=cell(num_vec,1);
% non_zero=cell(num_vec,1);
% powers_vector=cell(num_vec,1);
% for i=0:DIM-1
%     value{i*3+1}=[1,sigma,-sigma].';
%     value{i*3+2}=[1,1,-1,-1,epsilon].';
%     value{i*3+3}=[1,1,-beta].';
%     
%     non_zero(i*3+1)=3;
%     non_zero(i*3+2)=5;
%     non_zero(i*3+3)=3;
%     
%     powers_scalar{i*3+1}=[0,1,1
%         0,0,0];
%     powers_scalar{i*3+2}=[0,1,1,1,1; 0,1,0,0,0];
%     powers_scalar{i*3+3}=[0,1,1;0,0,0];
%     
%     powers_vector{i*3+1}=zeros(num_vec,non_zero{i*3+1});
%     powers_vector{i*3+1}(i*3+1,1)=1;
%     powers_vector{i*3+1}(i*3+2,2)=1;
%     powers_vector{i*3+1}(i*3+1,3)=1;
%     
%     powers_vector{i*3+2}=zeros(num_vec,non_zero{i*3+2});
%     powers_vector{i*3+2}(i*3+2,1)=1;
%     powers_vector{i*3+2}(i*3+1,2)=1;
%     powers_vector{i*3+2}(i*3+1,3)=1;
%     powers_vector{i*3+2}(i*3+3,3)=1;
%     %powers_vector{i*3+2}(i*3+1,5)=1;
%     %powers_vector{i*3+2}(i*3+2,4)=1;
%     
%     powers_vector{i*3+3}=zeros(num_vec,non_zero{i*3+3});
%     powers_vector{i*3+3}(i*3+3,1)=1;
%     powers_vector{i*3+3}(i*3+1,2)=1;
%     powers_vector{i*3+3}(i*3+2,2)=1;
%     powers_vector{i*3+3}(i*3+3,3)=1;
%     if i>0  % EPSILON term
%         powers_vector{i*3+2}((i-1)*3+1,5)=1;
%         powers_vector{i*3+2}((i-1)*3+2,4)=1;
%         
%         %powers_vector{i*3+1}((i-1)*3+1,3)=1;
%         %powers_vector{i*3+3}((i-1)*3+3,3)=1;
%     else
%         powers_vector{i*3+2}((DIM-1)*3+1,5)=1;
%         powers_vector{i*3+2}((DIM-1)*3+2,4)=1;
%         
%         %powers_vector{i*3+1}((DIM-1)*3+1,3)=1;
%         %powers_vector{i*3+3}((DIM-1)*3+3,3)=1;
%     end
%     
% end
% num_scal=2;
% deg_scal=2;deg_vec=2;
% alpha_coef=coefs(num_scal,num_vec,deg_scal,deg_vec,non_zero,...
%     powers_scalar,powers_vector,value);
% save(sprintf('%s_cont',s),'alpha_coef');
% 
% num_scal=1;
% for i=0:DIM-1
%     value{i*3+2}=[-2*pi*1i,pho,-1,-1,epsilon]/(2*pi);
%     powers_scalar{i*3+1}=[0,1,1];
%     powers_scalar{i*3+2}=[0,1,1,1,1];
%     powers_scalar{i*3+3}=[0,1,1];
% end
% num_scal=1;
% deg_scal=1;
% alpha_coef=coefs(num_scal,num_vec,deg_scal,deg_vec,non_zero,...
%     powers_scalar,powers_vector,value);
% save(s,'alpha_coef');
