function [s]=alpha_huge_lorenz(DIM, epsilon)
% function [s]=alpha_huge_lorenz(DIM, epsilon)
%
% INPUT
% DIM      number of correlated Lorenz systems
% epsilon  perturbation parameter (DEFAULT 10^-4) 
% OUTPUT
% s        name of the saved system
%
% alpha for super_huge Lorenz:
% DIM correlated Lorenz systems
% 
% (x1,y1,z1)...(xDIM,yDIM,zDIM)
%
% each solving \dot( x _i,y_i,z_i) =L + epsilon x_(i-1)
% considering x_0=x_7


s=sprintf('huge_lorenz%d',DIM);

num_vec=DIM*3;
sigma=10;beta=8/3;pho=28;
if nargin==1
    epsilon=10^-4;
end

value=cell(num_vec,1);
powers_scalar=cell(num_vec,1);
non_zero=cell(num_vec,1);
powers_vector=cell(num_vec,1);
for i=0:DIM-1
    value{i*3+1}=[-2*pi*1i,sigma,-sigma].'/(2*pi);
    value{i*3+2}=[-2*pi*1i,1,-1,-1,epsilon].'/(2*pi);
    value{i*3+3}=[-2*pi*1i,1,-beta].'/(2*pi);
    
    non_zero{i*3+1}=3;
    non_zero{i*3+2}=5;
    non_zero{i*3+3}=3;
    
    powers_scalar{i*3+1}=[0,1,1
        0,0,0];
    powers_scalar{i*3+2}=[0,1,1,1,1; 0,1,0,0,0];
    powers_scalar{i*3+3}=[0,1,1;0,0,0];
    
    powers_vector{i*3+1}=zeros(num_vec,non_zero{i*3+1});
    powers_vector{i*3+1}(i*3+1,1)=1;
    powers_vector{i*3+1}(i*3+2,2)=1;
    powers_vector{i*3+1}(i*3+1,3)=1;
    
    powers_vector{i*3+2}=zeros(num_vec,non_zero{i*3+2});
    powers_vector{i*3+2}(i*3+2,1)=1;
    powers_vector{i*3+2}(i*3+1,2)=1;
    powers_vector{i*3+2}(i*3+1,3)=1;
    powers_vector{i*3+2}(i*3+3,3)=1;
    %powers_vector{i*3+2}(i*3+1,5)=1;
    %powers_vector{i*3+2}(i*3+2,4)=1;
    
    powers_vector{i*3+3}=zeros(num_vec,non_zero{i*3+3});
    powers_vector{i*3+3}(i*3+3,1)=1;
    powers_vector{i*3+3}(i*3+1,2)=1;
    powers_vector{i*3+3}(i*3+2,2)=1;
    powers_vector{i*3+3}(i*3+3,3)=1;
    if i>0  % EPSILON term
        powers_vector{i*3+2}((i-1)*3+1,5)=1;
        powers_vector{i*3+2}((i-1)*3+2,4)=1;
        
        %powers_vector{i*3+1}((i-1)*3+1,3)=1;
        %powers_vector{i*3+3}((i-1)*3+3,3)=1;
    else
        powers_vector{i*3+2}((DIM-1)*3+1,5)=1;
        powers_vector{i*3+2}((DIM-1)*3+2,4)=1;
        
        %powers_vector{i*3+1}((DIM-1)*3+1,3)=1;
        %powers_vector{i*3+3}((DIM-1)*3+3,3)=1;
    end
    
end
num_scal=2;
deg_scal=2;deg_vec=2;
alpha_coef=coefs(num_scal,num_vec,deg_scal,deg_vec,non_zero,...
    powers_scalar,powers_vector,value);
save(sprintf('%s_cont',s),'alpha_coef');

num_scal=1;
for i=0:DIM-1
    value{i*3+2}=[-2*pi*1i,pho,-1,-1,epsilon]/(2*pi);
    powers_scalar{i*3+1}=[0,1,1];
    powers_scalar{i*3+2}=[0,1,1,1,1];
    powers_scalar{i*3+3}=[0,1,1];
end
num_scal=1;
deg_scal=1;
alpha_coef=coefs(num_scal,num_vec,deg_scal,deg_vec,non_zero,...
    powers_scalar,powers_vector,value);
save(s,'alpha_coef');