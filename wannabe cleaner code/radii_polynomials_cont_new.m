function [flag,Imin,Imax,new_iter,Yvector,Z0vector,Z1vector,Z2vector,new_step,Ys] = radii_polynomials_cont_new(xBar0,xBar1,DH0,DH1,...
    alpha0, alpha1,previous_iter,A0,A1)
% INPUT:
% xBar0      Xi_vector, approximate solution of the problem in x_0;
% xBar1      Xi_vector, approximate solution of the problem in x_1;
% alpha0,alpha1      full_problem, coefficients of the ODE equations;
% DH0        derivative in x_0;
% DH1        derivative in x_1;
%
% OUTPUT:
% flag      0  - program failed
%           1  - all good
%           2  - really good (consider decreasing delta)
% Imin      positive real value, left bound of the interval;
% Imax      positive real value, right bound of the interval;
% new_iter  storage of elements that are useful for next computation
% Yvector   value of the Y bound
% Z0vector  value of the Z0 bound
% Z1vector  value of the Z1 bound
% Z2vector  value of the Z2 bound
% new_step  euristic approximation of the new stepsize
%
% In this function, the radii polynomial refering to the given problem and
% solutions are computed.

global use_intlab
global nu
global Display
global talkative

% ERROR Spotting
if any(size(DH0)~=size(DH1))
    error('Sizes of derivatives must be equal.')
end
if length(size(DH0))>2
    error('Derivatives must be matrices.')
end
if size(DH0,1)~=size(DH0,2)
    error('Derivatives must be square matrices.');
end

if xBar0.size_scalar~=xBar1.size_scalar
    error('Scalar dimensions are not consistent.')
end
if xBar0.size_vector~=xBar1.size_vector
    error('Vector dimensions are not consistent.')
end
if  ~compatible_vec(alpha0,xBar0) || ~compatible_vec(alpha1,xBar0) || ~compatible_vec(alpha1,xBar1)
    error('Problems and solutions not compatible.')
end
if ~square(alpha1)
    error('Problem is not square, a unique solution does not exist.')
end
if nargin<6
    error('Too few input arguments')
end
if nargin<7
    previous_iter = [];
end


% default retunr in case of crash
flag=0;
Imin=[];
Imax=[];
Yvector=[];Z0vector=[];Z1vector=[];Z2vector=[];
new_iter=struct('Y',[],'Z1',[]);
new_step=norm(norm(xBar1-xBar0));

% important values
if ~exist('A0','var') || isempty(A0)
    A0=inv(DH0);
    if talkative>1
        fprintf('Computed A0, time %s\n',datestr(now,13));
    end
end
if any(size(A0)~=size(DH0))
    warning('A0 dimension not matching with DH0 dimension');
    A0=inv(DH0);
    if talkative>1
        fprintf('Computed A0, time %s\n',datestr(now,13));
    end
end


if ~exist('A1','var') || isempty(A1)
    A1=inv(DH1);
    if talkative>1
        fprintf('Computed A1, time %s\n',datestr(now,13));
    end
end
if any(size(A1)~=size(DH1))
    warning('A0 dimension not matching with DH0 dimension');
    A1=inv(DH1);
    if talkative>1
        fprintf('Computed A1, time %s\n',datestr(now,13));
    end
end


Adagger_delta=DH1-DH0;

% symmetrise A0 and A1
%  A0=symmetrise_A(A0,xBar0);
%  A1=symmetrise_A(A1,xBar0);

% change into intvals
if use_intlab
    DH0=intval(DH0);
    DH1=intval(DH1);
    A0=intval(A0);
    A1=intval(A1);
    Adagger_delta=intval(Adagger_delta);
    nu=intval(nu);
    xBar0=intval(xBar0);
    xBar1=intval(xBar1);
    for i=1:xBar0.size_vector
        alpha0.vector_field.value{i}=intval(alpha0.vector_field.value{i});
        alpha1.vector_field.value{i}=intval(alpha1.vector_field.value{i});
    end
    for i = 1:alpha0.scalar_equations.number_equations_pol
        alpha0.scalar_equations.polynomial_equations.value{i}=intval(alpha0.scalar_equations.polynomial_equations.value{i});
        alpha1.scalar_equations.polynomial_equations.value{i}=intval(alpha1.scalar_equations.polynomial_equations.value{i});
    end
    %for i=1:2
    %    coefs_linear{i}=intval(coefs_linear{i});
    %end
end

try
    if ~ isempty(previous_iter) && ~isempty(previous_iter.Y)
        [Yvector,new_iter_Y,Ys]=Y_bound_cont_new(A0,A1,xBar0,xBar1,alpha0,alpha1,previous_iter.Y);
    else
        [Yvector,new_iter_Y,Ys]=Y_bound_cont_new(A0,A1,xBar0,xBar1,alpha0,alpha1);
    end
    
    if talkative>1
        fprintf('\nComputed Y, time %s\n\n',datestr(now,13));
    end
catch
    return
end

try
    [Z0vector,Z0s]=Z0_bound_cont(DH0,DH1,A0,A1,xBar0,xBar1);
    
    if talkative>1
        fprintf('\nComputed Z0, time %s\n\n',datestr(now,13));
    end
catch
    return
end


try
    if  ~ isempty(previous_iter) && ~isempty(previous_iter.Z1)
        [Z1vector,new_iter_Z1,Z1s]=Z1_bound_cont_new(A0,A1,xBar0,xBar1,alpha0,...
            alpha1,Adagger_delta,previous_iter.Z1);
    else
        [Z1vector,new_iter_Z1,Z1s]=Z1_bound_cont_new(A0,A1,xBar0,xBar1,alpha0,...
            alpha1,Adagger_delta);
    end
    if talkative>1
        fprintf('\nComputed Z1, time %s\n\n',datestr(now,13));
    end
catch
    return
end

[Z2vector,Z2s]=Z2_bound_cont_new(A0,A1,xBar0,xBar1,alpha0);

    if talkative>1
        fprintf('\nComputed Z2, time %s\n\n',datestr(now,13));
    end
new_iter=struct('Y',new_iter_Y,'Z1',new_iter_Z1);

try
    [Imin,Imax]=find_negative(Z2vector,Z1vector,Z0vector,Yvector);
    if talkative>1
        fprintf('\n Computed validation interval, time %s\n\n',datestr(now,13));
    end
catch
    return
end

if any(size(Imin))<1 || min(Imin,Imax)<0
    return
end

% computation of ideal stepsize
r0=(Imin+Imax)/2;
poly_s=Ys+r0*(Z0s+Z1s)+r0^2*Z2s;
poly_s(:,end)=poly_s(:,end)-r0;
while sum(abs(poly_s(:,1)))==0
    poly_s=poly_s(:,2:end);
end
s_possible=max(roots(max(poly_s)));
%s_possible=zeros(xBar0.size_vector+xBar0.size_scalar,length(poly_s(1,:))-1);
%for i=1:xBar0.size_vector+xBar0.size_scalar
%    s_possible(i,:)=roots(poly_s(i,:));
%end
%s_possible=min(abs(s_possible(imag(s_possible)==0))-1);
if s_possible>.15
    new_step=(1+s_possible)/2;
elseif s_possible>.1
    new_step=1;
elseif s_possible>0.3
    new_step=0.95;
else
    new_step=0.95-abs(s_possible);
end

% plot
if Display
    p=@(r) Z0vector*r+Z1vector*r+Z2vector*r.^2-ones(size(Z0vector))*r+Yvector*ones(size(r));
    xplot=0:0.000001:Imax*1.5;
    figure;
    pplot=p(xplot);
    plot(xplot,pplot,xplot,0*xplot,'k');
    hold on
    plot(Imin,0,'k*',Imax,0,'k*')
end

if Imax-Imin>1e-4
    flag=2;
else
    flag=1;
end

return
end





function A_sim=symmetrise_A(A,x)
A_sim=A;

%A_sim(1:x.size_scalar,1:x.size_scalar)=(A_sim(1:x.size_scalar,1:x.size_scalar)+...
%    conj(A_sim(1:x.size_scalar,1:x.size_scalar)))/2;

length_nodes=2*x.nodes+1;
sym_vec=@(v) ( v + conj(v(end:-1:1,:)))/2;
sym_mat=@(A) ( A + conj(A(end:-1:1,end:-1:1)))/2;

for j=1:x.size_vector
    vec=x.size_scalar+(j-1)*length_nodes+1:x.size_scalar+j*length_nodes;
    all_scal=1:x.size_scalar;
    %A_sim(all_scal,vec)=sym_vec(A_sim(all_scal,vec));
    A_sim(vec,all_scal)=sym_vec(A_sim(vec,all_scal));
end

for jj=1:x.size_vector
    for kk=1:x.size_vector
        vec1=x.size_scalar+(jj-1)*length_nodes+1:x.size_scalar+jj*length_nodes;
        vec2=x.size_scalar+(kk-1)*length_nodes+1:x.size_scalar+kk*length_nodes;
        A_sim(vec1,vec2)=sym_mat(A_sim(vec1,vec2));
    end
end


end