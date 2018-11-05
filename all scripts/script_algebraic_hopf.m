% the code request you to specify fn_  derivatives_ second_der and
% third_der
% the codes has been completely generalized to N dimensions, no troubles
% with that
% the codes now run on the 4 dim, hypercahotic Rossler system (no periodic
% orbit continuation just yet)
% works as of 4th July 2108

%addpath('./Intlab_V7.1/');
%startintlab

% assumption: length(X) = 2 + 3 length(x)

%f=@fn_Brusselator; % We choose a map.
%df=@derivatives_Brusselator;
%ddf=@second_der_Brusselator;
%dddf=@third_der_Brusselator;
% N=2;
% DIM = N;
% phi=[0.814284826068816
%   0.243524968724989]+ 0.0001*rand(2,1);


%f=@fn_Bautin; % We choose a map.

% f = @fn_Rossler; % We choose a map.
% df = @derivatives_Rossler;
% ddf = @second_der_Rossler;
% dddf = @third_der_Rossler;
% N = 4;
% DIM= N; % works mostly (still random guess)



% ALL DATA FOR LORENZ 84
f = @fn_Lorenz84; % We choose a map.
df = @derivatives_Lorenz84;
ddf = @second_der_Lorenz84;
dddf = @third_der_Lorenz84;
N = 4;
DIM= N;
X = [0.02 0.52   0.005   0.6   0.91   0.2   0.8976   0.7768   0.72   0.155   0.538   0.2578   0.47   0.51]'; %lorenx84
phi = [0.347693611712440   0.362607042331847   0.999178390839916   0.649549106979193]'; %lorenx84


% f = @fn_hyper; % We choose a map.
% df = @derivatives_hyper;
% ddf = @second_der_hyper;
% dddf = @third_der_hyper;
% N = 4;
% DIM= N;
% phi = [0.905682106663319
%    0.381464175365932
%    0.664833832475170
%    0.368689083748346];
% X = [0.922263414809001
%    0.394970627320998
%    0.881757829928789
%    0.957770670583953
%    0.951563501679421
%    0.067613225413794
%    0.433091257190636
%    0.153359080393839
%    0.029920996267439
%    0.830344092041568
%    0.460020830546845
%    0.131807532006908
%    0.962522674607118
%    0.166445597032489];

% for a=0.1, x=[-0.060578090252706   0.002256024778798  -0.090978632950962
% -0.040896226254254] has a periodic orbit of period 4.3


        % X0=rand(2+3*N,1); 
        % X=X0;
% end  0.0568   -0.5300    1.1975   -0.0335    0.2032   -0.4003   -2.0312    0.1300    0.5707    0.1367    0.7407   -0.3338   -0.1440    1.5509

%phi=[0.814284826068816
%   0.243524968724989];%

        % phi = rand(N,1);%ones(4,1);%[0;1;1;0];%
[X]=newton_Hopf(f,X,phi);

alpha=X(1); % Parameter 
beta=X(2); % distance of the complex conjugate eigs from the real axis

x=X(2+(1:N)); % the variables from the model
v1=X(N+2+(1:N)); v2=X(2*N+2+(1:N)); % The real and imaginary parts of the eigenvector

%%%% validation
% for the validation, the exact derivative of f is necessary
DF = big_derivative(df,X,phi);
if rank(DF)-size(DF,1)<0
    error('Derivative not invertible')
end

A=inv(DF);

F=F_general_Hopf(f,X,phi,df);

% start intlab
Y = norm(A*F);

Z1 = norm( eye(length(X)) - A *DF);

R = 0.3; % first approximation of the maximum radius

%DDF = second_derivative_Brusselator(X,phi,R);

XandR=infsup(X-R,X+R);
DDF= second_derivative_F(XandR,phi,df);

Z2 = norm(sup(abs(A*DDF)));

pol=@(r) Y+(Z1-1)*r+Z2*r.^2;

delta= (intval(Z1)-1)^2-4*intval(Z2)*intval(Y);
if delta<=0 
    error('Hopf not validated')
end
rmin=sup((-(Z1-1)-sqrt(delta))/(2*Z2));
rmax=inf((-(Z1-1)+sqrt(delta))/(2*Z2));

% check rmax<R 

if rmax>R
    warning('computation has to be rerun with higher R')
end

r=linspace(0,1.3*rmax,100);
plot(r,pol(r),'b',rmin,0,'ok',rmax,0,'or');


% compute the first Lyapunov coefficient as in Kuznetsov
q =infsup( v1 + 1i*v2-rmin, v1 + 1i*v2+rmin);

% p: DF(xH)^Tp=conj(beta)p
% p = null( DF(xH)^T-conj(beta)I)

df_mat=df(x,alpha);
% next line does not work
%p = null( df_mat.'-conj(1i*beta)*eye(size(df_mat)));

[all_p,all_beta]=eigs(df_mat.');
all_eigs=diag(all_beta);
[~,index_beta]=min(all_eigs-conj(1i*beta));
p=all_p(:,index_beta);
% verification of the eigenvalue
[l,p]=verifyeig(midrad(df_mat',rmin),conj(1i*beta),p);

complex_product=@(a,b) sum(conj(a).*b);


% look for k, the rescaling coefficient such that <kp,q>=1
% k = k1 + i k2

% k1/k2
ratio = real(complex_product(p,q))/imag(complex_product(p,q));
k2 = 1/( ratio* real(complex_product(p,q))+imag(complex_product(p,q)));
k1= k2*ratio;

p=(k1+1i*k2)*p;

beta_ver=infsup(beta-rmin,beta+rmin);

% compute the first lyapunov coefficient
l1 = 1/(2*beta_ver) * real(1i*complex_product(p,ddf(zeros(N,1),q,q)*complex_product(p,ddf(zeros(N,1),q,conj(q))))+...
    beta_ver*complex_product(p,dddf(zeros(N,1),q,q,conj(q))));

% check the FLC is nonzero
if intval(inf(l1))*intval(sup(l1))<0
    error('Could not validate Hopf, since the first Lyapunov coefficient is not verifyed to be non-zero')
end


% last check: all ther eigenvalues different then zero
[all_eigenvectors,all_eigenvalues]=eigs(df_mat.');
all_eigenvalues = diag(all_eigenvalues);
[~,index_beta1]=min(all_eigenvalues-(1i*beta));
[~,index_beta2]=min(all_eigenvalues+(1i*beta));
for i=1:length(all_eigs)
    if i~=index_beta1 && i~=index_beta2
        [L,~] = verifyeig(df_mat,all_eigenvalues(i),all_eigenvectors(:,i));
        if intval(inf(real(L)))*intval(sup(real(L)))<0
            error('Could not validate Hopf, since there is at least one eigenvalue that not verifyed to be non-zero')
        end
    end

end
fprintf('SUCCESS\n The Hopf bifurcation at [%1.3f,%1.3f], %1.3f has been verified\n',x(1),x(2),alpha)
fprintf('with a radius of %e to %e\n\n',rmin, rmax) 

x_star = X(3:6);
lambda_star = X(1);
eigenvec = v1+1i*v2;
eigenval = 1i*beta;
save('hopf_in_Rossler', 'X','l1','x_star','lambda_star','eigenvec','eigenval');



return



