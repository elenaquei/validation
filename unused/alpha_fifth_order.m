delta_value=0.2;
mu=1;

value=cell(2,1);
value{1}=[-2*pi*1i,1,-1,mu,1].'/(2*pi);
value{2}=[-2*pi*1i,-1].'/(2*pi);
powers_scalar=cell(2,1);
powers_scalar{1}=[0,1,1,1,1;0,0,0,0,1];
powers_scalar{2}=[0,1;0,0];

powers_vector=cell(2,1);
powers_vector{1}=[1,0,5,3,1;0,1,0,0,0];
powers_vector{2}=[0,1;1,0];
non_zero=cell(2,1);
non_zero{1}=5;non_zero{2}=2;
alpha_coef=coefs(2,2,2,5,non_zero,powers_scalar,powers_vector,value);
save('rychkov_cont','alpha_coef');


value=cell(2,1);
value{1}=[-2*pi*1i,1,-1,mu,-delta_value].'/(2*pi);
value{2}=[-2*pi*1i,-1].'/(2*pi);
powers_scalar=cell(2,1);
powers_scalar{1}=[0,1,1,1,1];
powers_scalar{2}=[0,1];

powers_vector=cell(2,1);
powers_vector{1}=[1,0,5,3,1;0,1,0,0,0];
powers_vector{2}=[0,1;1,0];
alpha_coef=coefs(1,2,2,5,non_zero,powers_scalar,powers_vector,value);
save('rychkov','alpha_coef');