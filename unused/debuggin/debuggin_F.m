% % % % debugging F and DF
% % % % stupid test
% % % % that works correctly
% % % % 
% % % % x1_dot + lambda * x1 * x2^2
% % % % lambda * x2_dot + x1^2 * x2
% % % %
% % % % but Newton keeps not converging: therefore there must be another mistake
% % % % somewhere!!!! D:
% % % 
% % % 
% % % nodes = 4;
% % % lin_coef= cell(3,1);
% % % lin_coef{1} = 0;
% % % lin_coef{2} = ones(1,2,nodes*2+1);
% % % lin_coef{3} = 1;
% % % pol = polynomial_coefs(1,2,0,[],[],[],[],[],[]);
% % % scalar = scalar_eq(1,0,1,2,lin_coef,pol);
% % % 
% % % dim_vec = 2;
% % % dim_scal = 1;
% % % n_eq = dim_vec;
% % % non_zero=[2,2];
% % % power = cell(2,1);
% % % power{1}{2} = [1,2]';
% % % power{1}{1} =[0 0]';
% % % 
% % % power{2}{1} = [0 0]';
% % % power{2}{2} = [2 1]';
% % % 
% % % power_lambda = cell(2,1);
% % % power_lambda{1} = [0 1];
% % % power_lambda{2} = [1, 0];
% % % coef = cell(2,1);
% % % coef{1} = [1 1];
% % % coef{2} = [1,1];
% % % der = power;
% % % der{1}{2} = [0 0]';
% % % for i =1:2
% % %     der{2}{i} =0* der{1}{2};
% % % end
% % % der{1}{1} = [1 0]';
% % % der{2}{1} = [0 1]';
% % % 
% % % vector = polynomial_coefs(dim_scal, dim_vec, n_eq, non_zero, coef, power_lambda, power, der);
% % % F = full_problem(scalar, vector);
% % % 
% % % n_nodes = 4;
% % % 
% % % sin_four = zeros(1,2*n_nodes+1);
% % % cos_four = sin_four;
% % % cos_four(n_nodes) = 1/2; cos_four(n_nodes+2) = 1/2;
% % % sin_four(n_nodes) = 1i/2; sin_four(n_nodes+2) = -1i/2;
% % % 
% % % y = [sin_four; cos_four];
% % % lambda = 7;
% % % sol = Xi_vector(lambda,y);
% % % residual = apply(F, sol);
% % % 
% % % exact1 = 1i*(-nodes:nodes) .* sin_four +  lambda*conv(conv(sin_four,cos_four,'same'),cos_four,'same');
% % % exact2 =  1i*(-nodes:nodes) .*lambda.*cos_four + conv(conv(sin_four,sin_four,'same'),cos_four,'same');
% % % 
% % % sum(abs(exact1- residual.vector(1,:)));
% % % sum(abs(exact2- residual.vector(2,:)));
% % % 
% % % DF = derivative(F,sol);
% % % 
% % % one_four(1:4*nodes+1) = 0;
% % % one_four(2*n_nodes+1) = 1;
% % % 
% % % exact11 = lambda * conv(cos_four,cos_four);
% % % exact12 = 2 * lambda*conv(sin_four,cos_four);
% % % exact21 = 2*conv(sin_four,cos_four);
% % % exact22 = lambda * one_four*0 + conv(sin_four,sin_four);
% % % 
% % % DF_int = DF.derivative_Fx_toeplix;
% % % DF_int11 = squeeze(DF_int(1,1,n_nodes+1:end-n_nodes));
% % % max(abs(DF_int11 -exact11.'));
% % % DF_int12 = squeeze(DF_int(2,1,n_nodes+1:end-n_nodes));
% % % max(abs(DF_int12 -exact12.'));
% % % DF_int21 = squeeze(DF_int(1,2,n_nodes+1:end-n_nodes));
% % % max(abs(DF_int21 -exact21.'));
% % % DF_int22 = squeeze(DF_int(2,2,n_nodes+1:end-n_nodes));
% % % max(abs(DF_int22 - exact22.'));
% % % 
% % % exact_lambda1 = conv(sin_four,conv(cos_four, cos_four));
% % % exact_lambda2 = 1i*(-3*nodes:3*nodes).*conv(cos_four,one_four);
% % % max(abs(exact_lambda1.' - squeeze(DF.derivative_Flambda(1,1,:))));
% % max(abs(exact_lambda2.' - squeeze(DF.derivative_Flambda(1,2,:))));
% % 
% % zero_four(1:2*nodes+1) =0;
% % pad =[zero_four,zero_four(1:end-2)];
% % big_pad = zeros(2,(25-9)/2);
% % temp =(squeeze(F.scalar_equations.linear_coef{2}));
% % lin_coef = [F.scalar_equations.linear_coef{1} reshape([big_pad,temp,big_pad].',[],1).'];
% % 
% exact_mat = [lin_coef;
%     exact_lambda1.' diag(-3*nodes:3*nodes)+toeplitz([exact11(2*nodes+1:-1:1),pad],[exact11(2*nodes+1:end),pad])   toeplitz([exact12(2*nodes+1:-1:1),pad],[exact12(2*nodes+1:end),pad])
%     exact_lambda2.' toeplitz([exact21(2*nodes+1:-1:1),pad],[exact21(2*nodes+1:end),pad])   lambda*diag(-3*nodes:3*nodes)+toeplitz([exact22(2*nodes+1:-1:1),pad],[exact22(2*nodes+1:end),pad])];
% 
% DF_mat = derivative_to_matrix(DF, DF.nodes);
% 
% max(max(abs(exact_mat- DF_mat)));
% 


% debugging F and DF
% stupid test
% 
% x1 + lambda * x1 * x2^2
% lambda * x2 + x1^2 * x2
%
% and apply Taylor to Hopf

nodes = 4;
n_nodes = nodes;
lin_coef= cell(3,1);
lin_coef{1} = 0;
lin_coef{2} = ones(1,2,nodes*2+1);
lin_coef{3} = 1;
pol = polynomial_coefs(1,2,0,[],[],[],[],[],[]);
scalar = scalar_eq(1,0,1,2,lin_coef,pol);

dim_vec = 2;
dim_scal = 1;
n_eq = dim_vec;
non_zero=[2,2];
power = cell(2,1);
power{1}{2} = [1,2]';
power{1}{1} =[1 0]';

power{2}{1} = [0 1]';
power{2}{2} = [2 1]';

power_lambda = cell(2,1);
power_lambda{1} = [0 1];
power_lambda{2} = [1, 0];
coef = cell(2,1);
coef{1} = [1 1];
coef{2} = [1,1];
der = power;
der{1}{2} = [0 0]';
for i =1:2
    der{2}{i} =0* der{1}{2};
end
der{1}{1} = [1 0]';
der{2}{1} = [0 1]';

vector = polynomial_coefs(dim_scal, dim_vec, n_eq, non_zero, coef, power_lambda, power, der);
f = set_zero_derivative(vector);
F = Taylor_series_Hopf(f,n_nodes); % debugged



sin_four = zeros(1,2*n_nodes+1);
cos_four = sin_four;
cos_four(n_nodes) = 1/2; cos_four(n_nodes+2) = 1/2;
sin_four(n_nodes) = 1i/2; sin_four(n_nodes+2) = -1i/2;

y = [sin_four; cos_four];
x = [0.5, 0.75];
lambda = 7;
T=1; a = 100;
sol = Xi_vector([T, lambda, a, x],y);
residual = apply(F, sol);

y1 = sin_four; y2 = cos_four; 
x1 = x(1); x2 = x(2);
Conv = @(x,y) conv(x,y,'same');
der = @(x) 1i*(-nodes:nodes).*x;

exact1 = der(y1) - T*(lambda*x2^2*y1 + 2*lambda*x1*x2*y2 + 2*lambda*x2*a*Conv(y1,y2) +lambda*x1*a*Conv(y2,y2)+...
    +lambda*a*a*Conv(Conv(y1,y2),y2) +y1);
exact2 = der(y2) - T*(x1^2*y2 + 2*x1*x2*y1 + 2*x1*a*Conv(y1,y2) +x2*a*Conv(y1,y1)+...
    +a*a*Conv(Conv(y1,y2),y1) +lambda*y2);

sum(abs(exact1- residual.vector(1,:)));
sum(abs(exact2- residual.vector(2,:)));

exact_scal1 = x1 + lambda * x1 * x2^2;
exact_scal2 = lambda * x2 + x1^2 * x2;

abs(exact_scal1 - residual.scalar(end-1))
abs(exact_scal2 - residual.scalar(end))

DF = derivative(F,sol);

Df_exact = [x1*x2^2, 1+lambda*x2^2 ,2*lambda*x1*x2;
    x2, 2*x1*x2, lambda + x1^2];

abs(Df_exact - DF.derivative_Glambda(4:5,[2,4,5]))
% 
% one_four(1:4*nodes+1) = 0;
% one_four(2*n_nodes+1) = 1;
% 
% exact11 = lambda * conv(cos_four,cos_four);
% exact12 = 2 * lambda*conv(sin_four,cos_four);
% exact21 = 2*conv(sin_four,cos_four);
% exact22 = lambda * one_four*0 + conv(sin_four,sin_four);
% 
% DF_int = DF.derivative_Fx_toeplix;
% DF_int11 = squeeze(DF_int(1,1,n_nodes+1:end-n_nodes));
% max(abs(DF_int11 -exact11.'));
% DF_int12 = squeeze(DF_int(2,1,n_nodes+1:end-n_nodes));
% max(abs(DF_int12 -exact12.'));
% DF_int21 = squeeze(DF_int(1,2,n_nodes+1:end-n_nodes));
% max(abs(DF_int21 -exact21.'));
% DF_int22 = squeeze(DF_int(2,2,n_nodes+1:end-n_nodes));
% max(abs(DF_int22 - exact22.'));
% 
% exact_lambda1 = conv(sin_four,conv(cos_four, cos_four));
% exact_lambda2 = 1i*(-3*nodes:3*nodes).*conv(cos_four,one_four);
% max(abs(exact_lambda1.' - squeeze(DF.derivative_Flambda(1,1,:))));
% max(abs(exact_lambda2.' - squeeze(DF.derivative_Flambda(1,2,:))));
% 
% zero_four(1:2*nodes+1) =0;
% pad =[zero_four,zero_four(1:end-2)];
% big_pad = zeros(2,(25-9)/2);
% temp =(squeeze(F.scalar_equations.linear_coef{2}));
% lin_coef = [F.scalar_equations.linear_coef{1} reshape([big_pad,temp,big_pad].',[],1).'];
% 
% exact_mat = [lin_coef;
%     exact_lambda1.' diag(-3*nodes:3*nodes)+toeplitz([exact11(2*nodes+1:-1:1),pad],[exact11(2*nodes+1:end),pad])   toeplitz([exact12(2*nodes+1:-1:1),pad],[exact12(2*nodes+1:end),pad])
%     exact_lambda2.' toeplitz([exact21(2*nodes+1:-1:1),pad],[exact21(2*nodes+1:end),pad])   lambda*diag(-3*nodes:3*nodes)+toeplitz([exact22(2*nodes+1:-1:1),pad],[exact22(2*nodes+1:end),pad])];

DF_mat = derivative_to_matrix(DF, DF.nodes);

% max(max(abs(exact_mat- DF_mat)));

n_nodes = DF.nodes;
sin_four = zeros(1,2*n_nodes+1);
cos_four = sin_four;
cos_four(n_nodes) = 1/2; cos_four(n_nodes+2) = 1/2;
sin_four(n_nodes) = 1i/2; sin_four(n_nodes+2) = -1i/2;

y = [sin_four; cos_four];
x = [0.5, 0.75];
lambda = 7;
T=1; a = 100;
sol = Xi_vector([T, lambda, a, x],y);
Sol = Xi_vec2vec(sol);

fin_dif_der = 0*DF_mat;
h = 10^-6;
fx = Xi_vec2vec(apply(F,sol));
e = eye(length(Sol));
for i = 1:length(Sol)
    x_h = vec2Xi_vec(Sol+e(:,i)*h,sol);
    fin_dif_der(:,i) = (Xi_vec2vec(apply(F,x_h))-fx)/h;
end

one_four(1:2*n_nodes+1) = 0;
one_four(n_nodes+1) = 1;
y1 = sin_four; y2 = cos_four; 

D1F1 = - lambda*x2^2*T*one_four - 2*lambda*T*x2*a*y2 -lambda*T*a^2*Conv(y2,y2) - T*one_four;
max(abs(D1F1.' - squeeze(DF.derivative_Fx_toeplix(1,1,:))))
D2F1 = - 2*lambda *T*x1*x2*one_four - 2*lambda*T*x2*a*y1 - 2*lambda*T*x1*a*y2 - 2*lambda*T*a^2*Conv(y1,y2);
max(abs(D2F1.' - squeeze(DF.derivative_Fx_toeplix(2,1,:))))

%D1F2 = - 2*T*x1*x2 *one_four - 2*lambda*T*a*x1*y2 - 2*T*a*x2*y1 - 2*T*a^2*Conv(y1,y2);
D1F2 = - T*( 2*x1*x2*one_four + 2*x1*a*Conv(one_four,y2) +2*x2*a*Conv(one_four,y1)+...
    +a*a*2*Conv(y1,y2));
max(abs(D1F2.' - squeeze(DF.derivative_Fx_toeplix(1,2,:))))
D2F2 =- T*(x1^2*one_four + 2*x1*a*Conv(y1,one_four) +...
    +a*a*Conv(Conv(y1,one_four),y1) +lambda*one_four);
%-T*lambda*one_four-T*x1^2*one_four-2*lambda*T*a*x1*y1-T*a^2*Conv(y1,y1);
max(abs(D2F2.' - squeeze(DF.derivative_Fx_toeplix(2,2,:))))

Dlambda1 = zeros(2*n_nodes+1,5);
Dlambda2 = Dlambda1;

zero_four(1:2*nodes+1) =0;
pad =zeros(1,n_nodes);
big_pad = zeros(2,(25-9)/2);
temp =(squeeze(F.scalar_equations.linear_coef{2}));
%lin_coef = [F.scalar_equations.linear_coef{1} reshape([big_pad,temp,big_pad].',[],1).'];
lin_coef = zeros(5,5+(n_nodes*2+1)*2);
exact_mat = [lin_coef;
    Dlambda1 diag(-n_nodes:n_nodes)+toeplitz([D1F1(n_nodes+1:-1:1),pad],[D1F1(n_nodes+1:end),pad])   toeplitz([D2F1(n_nodes+1:-1:1),pad],[D2F1(n_nodes+1:end),pad])
    Dlambda2 toeplitz([D1F2(n_nodes+1:-1:1),pad],[D1F2(n_nodes+1:end),pad])   diag(-n_nodes:n_nodes)+toeplitz([D2F2(n_nodes+1:-1:1),pad],[D2F2(n_nodes+1:end),pad])];
