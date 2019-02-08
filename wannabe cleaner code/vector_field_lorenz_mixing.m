function [alpha_cont,s]=vector_field_lorenz_mixing(DIM,s)
% function [alpha_cont,s]=vector_field_lorenz_mixing(DIM,s)
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
% with intertwined dependency of the type
% x_i dot = sigma( y_i - x_i)
% y_i dot = x_i ( rho + eps*i - z_i ) - y_(i+1)
% z_i dot = x_i y _i - beta z_i
% where rho is constant and we continue in eps as long as possible.

if nargin<3 || isempty(s)
    s=sprintf('mixing_lorenz%d',DIM);
end

sigma=10;beta=8/3; pho=28;

string_lorenz = '- dot x_3i+1 + sigma l1 x_3i+2 - sigma l1 x_3i+1 \n - dot x_3i+2 + pho l1 x_3i+1 + I l1 l2 x_3i+1 - l1 x_3i+1 x_3i+3 - l1 x_3i+5 \n - dot x_3i+3 + l1 x_3i+1 x_3i+2 - beta l1 x_3i+3\n'; % general lorenz
string_lorenz_vars = strrep(string_lorenz, 'sigma' , num2str(sigma)); % plugging in sigma
string_lorenz_vars = strrep(string_lorenz_vars, 'beta' , num2str(beta)); % plugging in beta
string_lorenz_cont = strrep(string_lorenz_vars, 'pho', num2str(pho)); % for point wise first system, plugging in pho

string_vf='';
for i=1:DIM
    string_lorenz_i=strrep(string_lorenz_cont, '_3i+1', num2str(mod(3*(i-1)+1-1,3*DIM)+1));
    string_lorenz_i=strrep(string_lorenz_i, '_3i+2', num2str(mod(3*(i-1)+2-1,3*DIM)+1));
    string_lorenz_i=strrep(string_lorenz_i, '_3i+3', num2str(mod(3*(i-1)+3-1,3*DIM)+1));
    string_lorenz_i=strrep(string_lorenz_i, '_3i+5', num2str(mod(3*(i-1)+5-1,3*DIM)+1));
    string_lorenz_i=strrep(string_lorenz_i, 'I', num2str(i));
    string_vf=strcat(string_vf,string_lorenz_i);
end

alpha_cont = from_string_to_polynomial_coef(string_vf);

return