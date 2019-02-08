function [alpha_cont,s]=coupling_VDP(DIM,s)
% function [alpha_cont,s]=coupling_VDP(DIM,s)
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
% with intertwined dependency for van der Pol model:
% x dot = y_1
% y_i dot = mu (1-x^2)*(y_i+y_(i+1))/2-x


% working as of 9th July 2019

if nargin<3 || isempty(s)
    s=sprintf('vanderPol_mixing%d',DIM);
end

mu = 1;

string_vd1 = '- dot x1 + l1 x2 \n';
string_vdp2 = ' - dot x_i + half_mu l1 x_i + half_mu l1 x_i+1 - half_mu x1^2 x_i - half_mu l1 x1^2 x_i+1 - x1\n';
string_vdp_vars = strrep(string_vdp2, 'half_mu' , num2str(mu/2)); % plugging in mu/2

string_vf=string_vd1;
for i=2:DIM+1
    if i+1<=DIM
        string_vdp_i=strrep(string_vdp_vars, '_i+1', num2str(i+1));
    else
        string_vdp_i=strrep(string_vdp_vars, '_i+1', num2str(2));
    end
    string_vdp_i=strrep(string_vdp_i, '_i', num2str(i));
    string_vf=strcat(string_vf,string_vdp_i);
end

alpha_cont = from_string_to_polynomial_coef(string_vf);

return