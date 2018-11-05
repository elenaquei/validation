% alpha of stupid mixing van der Pol
%
% \dot x_1 
% \dot y_1 = V1
%
% \dot x_2
% \dot y_2 = V2 ( 1 + epsilon x1 y1)
%
% \dot x_3
% \dot y_3 = V3 ( 1 + epsilon x1^2 y1^2 x2^2 y2^2 )

epsilon = 0;
num_vec=4;

value=cell(num_vec,1);
value{1}=[-2*pi*1i,1].'/(2*pi);
value{2}=[-2*pi*1i,1, -1, -1].'/(2*pi);
value{3}=[-2*pi*1i,1,epsilon].'/(2*pi);
value{4}=[-2*pi*1i,1,-1,-1,epsilon*1,epsilon*1,epsilon].'/(2*pi);
%value{5}=[-2*pi*1i,1,epsilon].'/(2*pi);
%value{6}=[-2*pi*1i,1, -1, -1, epsilon*1, -epsilon*1, -epsilon*1].'/(2*pi);

powers_scalar=cell(num_vec,1);
powers_scalar{1}=[0,1;
    0,0];
powers_scalar{2}=[0,1,1,1;
    0,1,1,0];
powers_scalar{3}=[0,1,1;
    0,0,0];
powers_scalar{4}=[0,1,1,1,1,1,1;
    0,1,1,0,1,1,0];
%powers_scalar{5}=[0,1,1;
%    0,0,0];
%powers_scalar{6}=[0,1,1,1,1,1,1;
%    0,1,1,0,1,1,1];

non_zero=cell(num_vec,1);
non_zero{1}=2;
non_zero{2}=4;
non_zero{3}=3;
non_zero{4}=7;
%non_zero{5}=3;
%non_zero{6}=7;


powers_vector=cell(num_vec,1);
for i=1:num_vec
    powers_vector{i}=zeros(num_vec,non_zero{i});
end
powers_vector{1}(1:2,:)=[1,0;0,1];
powers_vector{2}(1:2,:)=[0,0,2,1;
    1,1,1,0];
powers_vector{3}(1:4,:)=[0,0,1
    0,0,1
    1,0,0
    0,1,1];
powers_vector{4}(1:4,:)=[0,0,0,0,1,1,1
    0,0,0,0,1,1,1
    0,0,2,1,0,2,1
    1,1,1,0,1,1,0];
%powers_vector{5}=[0,0,2
%    0,0,2
%    0,0,2
%    0,0,2
%    1,0,0
%    0,1,1];
%powers_vector{6}=[0,0,0,0,2,2,2
%    0,0,0,0,2,2,2
%    0,0,0,0,2,2,2
%    0,0,0,0,2,2,2
%    0,0,2,1,0,2,1
%    1,1,1,0,1,1,0];
alpha_coef=coefs(2,num_vec,2,9,non_zero,powers_scalar,powers_vector,value);
save('mixingVDP_cont','alpha_coef');

mu=1.2;
value=cell(num_vec,1);
value{1}=[-2*pi*1i,1].'/(2*pi);
value{2}=[-2*pi*1i,mu, -mu, -1].'/(2*pi);
value{3}=[-2*pi*1i,1,epsilon].'/(2*pi);
value{4}=[-2*pi*1i,mu,-mu,-1,epsilon*mu,epsilon*mu,epsilon].'/(2*pi);
value{5}=[-2*pi*1i,1,epsilon].'/(2*pi);
value{6}=[-2*pi*1i,mu, -mu, -1, epsilon*mu, -epsilon*mu, -epsilon*mu].'/(2*pi);

%powers_scalar=cell(num_vec,1);
for i=1:num_vec
    powers_scalar{i}=powers_scalar{i}(1,:);
end

alpha_coef=coefs(1,num_vec,2,9,non_zero,powers_scalar,powers_vector,value);
save('mixingVDP','alpha_coef');
