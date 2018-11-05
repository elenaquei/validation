function [x,x_dot, DH,xShort,alpha_coef,coefs_linear,alpha0,coefs_linearshort] = compute_solution_and_derivative(Name_system,flag_plot,n_nodes,maxiter,min_res)
% function [x,x_dot, DH,xShort,alpha_coef,coefs_linear] = 
%   compute_solution_and_derivative(Name_system,flag_plot,n_nodes,...
%   maxiter,min_res)
% 
% INPUT 
% Name_system    string that determines the system to be solved
% flag_plot      1 if you want plots, 0 of not
% n_nodes        requested number of nodes of the solution
% maxiter        maximum number of iterations for Newton
% min_res        requested residual for Newton
% 
% OUTPUT
% x             Xi_vector, solution of the full system (all variables included)
% x_dot         complex vector in the null-space of the non-square derivative
% DH            square matrix, full derivative
% xShort        Xi_vector, solution before the addition of the last variable
% alha_coefs    coefs, coefficients of the vector field
% coefs_linear  cell, coefficients for the scalar equations



% set the necessary variables inside a switch
switch Name_system
    case 'test_case'
        mu=1;
        new_var=mu;
        
        init_coord=[1,0];
        approx_period=2*pi;
        
        % FOURIER COEFS
        disp('sin-cos coefficients');
        Name_system='test_case.mat';
        % interactive_constructor is a script that recalls the already
        % stored coefficients of the vector field or interactively
        % construct a new vector field
        interactive_constructor
        
        alpha0=alpha_coef;
        
        coefs_linear0=coefs_linear;
        
        disp('')
        disp('continuation coefficients')
        Name_system='test_case_cont';
        % here te same function is called for the continuation system
        % (already including the powers of the parameter)
        interactive_constructor
        
        coefs_linear=coefs_linear0;
        
        alpha0.value{2}=[1i,-mu/(2*pi)];
        % this is due to the dependence
        % on mu, that needs to be computed on the run
        xBarGen0=solve_system(alpha0,n_nodes,init_coord,approx_period);
        % numerically solves the problem
        
        coefs_linear{3}(1)=0;
        for i=1:alpha0.size_scalar
            coefs_linear{3}(i)=-coefs_linear{1}(i,:)*xBarGen0.scalar.'-...
                sum(sum(squeeze(coefs_linear{2}(i,:,:)).*xBarGen0.vector,1));
        end
    case 'test_lorenz'
        sigma=10;
        beta=8/3;
        col=2;
        load LorNumericalSolutions
        pho=LorNumericalSolutions(1,col);
        
        non_zero{1}=3;
        non_zero{2}=4;
        non_zero{3}=3;
        
        value{1} = [1i, -sigma/(2*pi), sigma/(2*pi)];
        value{2} = [1i, -pho/(2*pi), 1/(2*pi), 1/(2*pi)];
        value{3} = [1i, -1/(2*pi), beta/(2*pi)];
        
        powers_scal{1}=[0, 1, 1];
        powers_scal{2}=[0,1,1,1];
        powers_scal{3}=[0,1,1];
        
        powers_vec{1}=[1,0,1;
            0,1,0;
            0,0,0];
        powers_vec{2}=[0,1,1,0;
            1,0,0,1;
            0,0,1,0];
        powers_vec{3}=[0,1,0;
            0,1,0;
            1,0,1];
        
        size_scal=1;
        size_vec=3;
        degscal=1;
        degvec=2;
        
        alpha0=coefs(size_scal,size_vec,degscal,degvec,...
            non_zero,powers_scal,powers_vec,value);
        
        
        
        value{2} = [1i, -1/(2*pi), 1/(2*pi), 1/(2*pi)];
        
        powers_scal{1}=[0, 1, 1;
            0,0,0];
        powers_scal{2}=[0,1,1,1;
            0,1,0,0];
        powers_scal{3}=[0,1,1;
            0,0,0];
        
        size_scal=2;
        degscal=2;
        
        alpha_coef=coefs(size_scal,size_vec,degscal,degvec,...
            non_zero,powers_scal,powers_vec,value);
        
        
        
        new_var = pho;
        x0=LorNumericalSolutions(3:5,col);
        x=zeros(3,n_nodes);
        
        for i=1:min((size(LorNumericalSolutions,1)-6/3),n_nodes)
            ck=LorNumericalSolutions(6+(i-1)*6:5+(i)*6,col);
            a=[ck(1),ck(3),ck(5)].';
            b=[ck(2),ck(4),ck(6)].';
            x(:,i)=a+1i*b;
        end
        for i=min((size(LorNumericalSolutions,1)-6/3),n_nodes):n_nodes
            x(:,i)=0;
        end
        
        x1=[conj(flip(x(1,1:end))),x0(1),x(1,:)];
        x2=[conj(flip(x(2,1:end))),x0(2),x(2,:)];
        x3=[conj(flip(x(3,1:end))),x0(3),x(3,:)];
        
        x=[x1;x2;x3];
        xBarGen0=Xi_vector(LorNumericalSolutions(2,col)/(2*pi),x);
        xBarGen0.scalar=xBarGen0.scalar^-1;
        
        coefs_linear0=cell(3,1);
        coefs_linear0{1}=zeros(alpha0.size_scalar,1);
        coefs_linear0{2}=zeros(alpha0.size_scalar,alpha0.size_vector,2*n_nodes+1);
        k0=0;
        coefs_linear0{2}(1,1,1+k0:end-k0)=1; % MOST SIMPLE CASE
        
        coefs_linear0{3}(1)=0;
        for i=1:alpha0.size_scalar
            coefs_linear0{3}(i)=-coefs_linear0{1}(i,:)*xBarGen0.scalar.'-...
                sum(sum(squeeze(coefs_linear0{2}(i,:,:)).*xBarGen0.vector,1));
        end
        coefs_linear=cell(3,1);
        coefs_linear{1}=zeros(alpha0.size_scalar,1);
        coefs_linear{2}=zeros(alpha0.size_scalar,alpha0.size_vector,2*n_nodes+1);
        k0=0;
        coefs_linear{2}(1,1,1+k0:end-k0)=1; % MOST SIMPLE CASE
        
        coefs_linear{3}(1)=0;
        for i=1:alpha0.size_scalar
            coefs_linear{3}(i)=-coefs_linear{1}(i,:)*xBarGen0.scalar.'-...
                sum(sum(squeeze(coefs_linear{2}(i,:,:)).*xBarGen0.vector,1));
        end
    case 'lorenz_cont'
        % LORENZ
        % clear all
        
        sigma=10;
        beta=8/3;
        col=2;
        load LorNumericalSolutions
        rho=LorNumericalSolutions(1,col);
        new_var = rho; % the variable we want to do continuation in
        x0=LorNumericalSolutions(3:5,col);
        x=zeros(3,n_nodes);
        
        for i=1:min((size(LorNumericalSolutions,1)-6/3),n_nodes)
            ck=LorNumericalSolutions(6+(i-1)*6:5+(i)*6,col);
            a=[ck(1),ck(3),ck(5)].';
            b=[ck(2),ck(4),ck(6)].';
            x(:,i)=a+1i*b;
        end
        for i=min((size(LorNumericalSolutions,1)-6/3),n_nodes):n_nodes
            x(:,i)=0;
        end
        
        x1=[conj(flip(x(1,1:end))),x0(1),x(1,:)];
        x2=[conj(flip(x(2,1:end))),x0(2),x(2,:)];
        x3=[conj(flip(x(3,1:end))),x0(3),x(3,:)];
        
        x=[x1;x2;x3];
        xBarGen0=Xi_vector(LorNumericalSolutions(2,col)/(2*pi),x);
        xBarGen0.scalar=xBarGen0.scalar^-1;
        
        disp('lorenz coefficients');
        Name_system='lorentz';
        % interactive_constructor is a script that recalls the already
        % stored coefficients of the vector field or interactively
        % construct a new vector field
        interactive_constructor
        
        alpha0=alpha_coef;
        
        coefs_linear0=coefs_linear;
        
        disp('')
        disp('continuation coefficients')
        Name_system='lorenz_cont';
        % here te same function is called for the continuation system
        % (already including the powers of the parameter)
        interactive_constructor
        
        coefs_linear=coefs_linear0;
        
        
        alpha0.value{1}=[-2*pi*1i,sigma,-sigma].'/(2*pi);
        alpha0.value{2}=[-2*pi*1i,rho,-1,-1].'/(2*pi);
        alpha0.value{3}=[-2*pi*1i,1,-beta].'/(2*pi);
        
        alpha_coef.value{1}=[-2*pi*1i,sigma,-sigma].'/(2*pi);
        alpha_coef.value{2}=[-2*pi*1i,1,-1,-1].'/(2*pi);
        alpha_coef.value{3}=[-2*pi*1i,1,-beta].'/(2*pi);
        % this is due to the dependence
        % on sigma, rho, beta, that needs to be computed on the run
        
        coefs_linear{3}(1)=0;
        for i=1:alpha0.size_scalar
            coefs_linear{3}(i)=-coefs_linear{1}(i,:)*xBarGen0.scalar.'-...
                sum(sum(squeeze(coefs_linear{2}(i,:,:)).*xBarGen0.vector,1));
        end
    case 'lorenz_christian_longer'
        disp('lorenz coefficients');
        
        sigma=10;
        beta=8/3;
        rho=28;
        new_var=rho;
        
        %init_coord=[-11.586951663722,-14.814615288122,27];
        %approx_period=3.7256417715558;
        init_coord=[-12.698941349915,-17.197497247713,27];
        approx_period=3.8695391125646;
        
        Name_system='lorentz';
        % interactive_constructor is a script that recalls the already
        % stored coefficients of the vector field or interactively
        % construct a new vector field
        interactive_constructor
        
        alpha0=alpha_coef;
        
        
        
        alpha0.value{1}=[-2*pi*1i,sigma,-sigma].'/(2*pi);
        alpha0.value{2}=[-2*pi*1i,rho,-1,-1].'/(2*pi);
        alpha0.value{3}=[-2*pi*1i,1,-beta].'/(2*pi);
        
        coefs_linear0=coefs_linear;
        
        xBarGen0=solve_system(alpha0,n_nodes,init_coord,approx_period);
        
        disp('')
        disp('continuation coefficients')
        Name_system='lorenz_cont';
        % here te same function is called for the continuation system
        % (already including the powers of the parameter)
        interactive_constructor
        
        coefs_linear=coefs_linear0;
        
        alpha_coef.value{1}=[-2*pi*1i,sigma,-sigma].'/(2*pi);
        alpha_coef.value{2}=[-2*pi*1i,1,-1,-1].'/(2*pi);
        alpha_coef.value{3}=[-2*pi*1i,1,-beta].'/(2*pi);
        % this is due to the dependence
        % on sigma, rho, beta, that needs to be computed on the run
        
        coefs_linear{3}(1)=0;
        for i=1:alpha0.size_scalar
            coefs_linear{3}(i)=-coefs_linear{1}(i,:)*xBarGen0.scalar.'-...
                sum(sum(squeeze(coefs_linear{2}(i,:,:)).*xBarGen0.vector,1));
        end
        
    case 'lorenz_christian'
        disp('lorenz coefficients');
        
        sigma=10;
        beta=8/3;
        rho=28;
        new_var=rho;
        
        init_coord=[-11.9985,-15.6842,27];
        approx_period=3.0235;
        
        Name_system='lorentz';
        % interactive_constructor is a script that recalls the already
        % stored coefficients of the vector field or interactively
        % construct a new vector field
        interactive_constructor
        
        alpha0=alpha_coef;
        
        
        
        alpha0.value{1}=[-2*pi*1i,sigma,-sigma].'/(2*pi);
        alpha0.value{2}=[-2*pi*1i,rho,-1,-1].'/(2*pi);
        alpha0.value{3}=[-2*pi*1i,1,-beta].'/(2*pi);
        
        coefs_linear0=coefs_linear;
        
        xBarGen0=solve_system(alpha0,n_nodes,init_coord,approx_period);
        
        disp('')
        disp('continuation coefficients')
        Name_system='lorenz_cont';
        % here te same function is called for the continuation system
        % (already including the powers of the parameter)
        interactive_constructor
        
        coefs_linear=coefs_linear0;
        
        alpha_coef.value{1}=[-2*pi*1i,sigma,-sigma].'/(2*pi);
        alpha_coef.value{2}=[-2*pi*1i,1,-1,-1].'/(2*pi);
        alpha_coef.value{3}=[-2*pi*1i,1,-beta].'/(2*pi);
        % this is due to the dependence
        % on sigma, rho, beta, that needs to be computed on the run
        
        coefs_linear{3}(1)=0;
        for i=1:alpha0.size_scalar
            coefs_linear{3}(i)=-coefs_linear{1}(i,:)*xBarGen0.scalar.'-...
                sum(sum(squeeze(coefs_linear{2}(i,:,:)).*xBarGen0.vector,1));
        end
        
    case 'lorenz_long'
        sigma=10;
        beta=8/3;
        load Orbit1.mat
        new_var=28;
        rho=28;
        disp('lorenz coefficients');
        Name_system='lorentz';
        % interactive_constructor is a script that recalls the already
        % stored coefficients of the vector field or interactively
        % construct a new vector field
        interactive_constructor
        
        alpha0=alpha_coef;
        
        coefs_linear0=coefs_linear;
        
        disp('')
        disp('continuation coefficients')
        Name_system='lorenz_cont';
        % here te same function is called for the continuation system
        % (already including the powers of the parameter)
        interactive_constructor
        coefs_linear=coefs_linear0;
        
        alpha_coef.value{1}=[-2*pi*1i,sigma,-sigma].'/(2*pi);
        alpha_coef.value{2}=[-2*pi*1i,1,-1,-1].'/(2*pi);
        alpha_coef.value{3}=[-2*pi*1i,1,-beta].'/(2*pi);
        % this is due to the dependence
        % on sigma, rho, beta, that needs to be computed on the run
        
        time=L;
        %y=y.';
        % yy2=y2(t2>=t2(end)-time,:);
        % y2=yy2;
        % t2=t2(t2>=t2(end)-time);
        % t2=t2-t2(1);
        
        %y=spline(t(1):0.0000001:t(end),y);
        y=spline(t,y,t(1):0.0000001:t(end));
        x=0*y;
        size_vec=3;
        for i=1:size_vec
            x(i,:)=1/(size(y,2))*fft(y(i,:));
        end
        xBar=cell(size_vec+alpha_coef.size_scalar,1);
        xBar{1}=time;
        m=n_nodes;
        for i=1:size_vec
            xBar{i+alpha_coef.size_scalar}=fftshift([x(i,1:m+1),x(i,end-m+1:end)]);
            xBar{i+alpha_coef.size_scalar}(1:m)=conj(flip(xBar{i+alpha_coef.size_scalar}(m+2:end)));
            xBar{i+alpha_coef.size_scalar}(m+1)=real(xBar{i+alpha_coef.size_scalar}(m+1));
        end
        
        size_scal=1;%alpha.size_scalar;
        nodes=m;
        %scal=xBar{1:alpha_coef.size_scalar};
        scal=[xBar{1}];%,rho];
        vec=[xBar{1+alpha_coef.size_scalar}];
        for i=2:size_vec
            vec=[vec;xBar{i+alpha_coef.size_scalar}];
        end
        xBarGen0=Xi_vector(scal,vec,1,size_vec,nodes);
        
        coefs_linear{3}(1)=0;
        for i=1:alpha0.size_scalar
            coefs_linear{3}(i)=... %-coefs_linear{1}(i,:)*xBarGen0.scalar.'-...
                sum(sum(squeeze(coefs_linear{2}(i,:,:)).*xBarGen0.vector,1));
        end
    case 'fifth_order'
        disp('Rychkov coefficients');
        
        mu=1;
        delta= 0.2;
        new_var=delta;
        
        init_coord=[0.893787,0];
        approx_period=2*pi;
        
        Name_system='rychkov';
        % interactive_constructor is a script that recalls the already
        % stored coefficients of the vector field or interactively
        % construct a new vector field
        interactive_constructor
        
        alpha0=alpha_coef;
        
        
        
        alpha0.value{1}=[-2*pi*1i,1,-1,mu,-delta].'/(2*pi);
        alpha0.value{2}=-[2*pi*1i,+1].'/(2*pi);
        
        coefs_linear0=coefs_linear;
        
        xBarGen0=solve_system(alpha0,n_nodes,init_coord,approx_period);
        
        disp('')
        disp('continuation coefficients')
        Name_system='rychkov_cont';
        % here te same function is called for the continuation system
        % (already including the powers of the parameter)
        interactive_constructor
        
        coefs_linear=coefs_linear0;
        
        alpha_coef.value{1}=[-2*pi*1i,1,-1,1,-1].'/(2*pi);
        alpha_coef.value{2}=-[2*pi*1i,1].'/(2*pi);
        % this is due to the dependence
        % on sigma, rho, beta, that needs to be computed on the run
        
        coefs_linear{3}(1)=0;
        for i=1:alpha0.size_scalar
            coefs_linear{3}(i)=-coefs_linear{1}(i,:)*xBarGen0.scalar.'-...
                sum(sum(squeeze(coefs_linear{2}(i,:,:)).*xBarGen0.vector,1));
        end
        
    case {'mixingVDP_cont','stupidmixing3VDP_cont'}
        % compute van der pol solution and triplicate it
        this_case=1;
        if strcmp(Name_system,'stupidmixing3VDP_cont')
            this_case=2;
        end
        mu=1.2;
        new_var=mu;
        
        alpha_mixingVDP;
        
        init_coord=[-1.288,0.9345];
        approx_period=6.3;
        
        disp('van der pol coefficients');
        Name_system='vanderpol';
        interactive_constructor
        
        alpha0=alpha_coef;
        
        xBarGen0_small=solve_system(alpha0,n_nodes,init_coord,approx_period);
        [xBarGen0_small]=Newton_Xi(xBarGen0_small,alpha0,coefs_linear,maxiter,min_res);
        
        disp('full mixing van der pol coefficients')
        if this_case==1
        Name_system= 'mixing3VDP';
        elseif this_case==2
            Name_system='stupidmixing3VDP';
        else 
            error('somthing wrong, why here?')
        end
        interactive_constructor;
        
        
        vector_big=zeros(alpha_coef.size_vector,2*xBarGen0_small.nodes+1);
        for i=0:alpha_coef.size_vector/2-1
            vector_big(i*2+(1:2),:)=xBarGen0_small.vector;
        end
        
        xBarGen0=Xi_vector(xBarGen0_small.scalar, vector_big,1,...
            alpha_coef.size_vector,xBarGen0_small.nodes);
        
        alpha0=alpha_coef;
        coefs_linear0=coefs_linear;
        
        disp('')
        disp('continuation coefficients')
        if this_case==1
            Name_system= 'mixing3VDP_cont';
        elseif this_case==2
            Name_system='stupidmixing3VDP_cont';
        else 
            error('somthing wrong, why here?')
        end
        % here te same function is called for the continuation system
        % (already including the powers of the parameter)
        interactive_constructor
        
        coefs_linear=coefs_linear0;
        
        coefs_linear{3}(1)=0;
        for i=1:alpha0.size_scalar
            coefs_linear{3}(i)=-coefs_linear{1}(i,:)*xBarGen0.scalar.'-...
                sum(sum(squeeze(coefs_linear{2}(i,:,:)).*xBarGen0.vector,1));
        end
        
    case 'vanderpol_cont'
        mu=1.2;
        new_var=mu;
        
        init_coord=[-1.288,0.9345];
        approx_period=6.3;
        
        % FOURIER COEFS
        disp('van der pol coefficients');
        Name_system='vanderpol';
        % interactive_constructor is a script that recalls the already
        % stored coefficients of the vector field or interactively
        % construct a new vector field
        interactive_constructor
        
        alpha0=alpha_coef;
        %alpha1=alpha_coef;
        coefs_linear0=coefs_linear;
        
        disp('')
        disp('continuation coefficients')
        Name_system='vanderpol_cont';
        % here te same function is called for the continuation system
        % (already including the powers of the parameter)
        interactive_constructor
        
        coefs_linear=coefs_linear0;
        
        alpha0.value{2}=[-2*pi*1i,mu,-mu,-1]/(2*pi); 
        % this is due to the dependence
        % on mu, that needs to be computed on the run
        xBarGen0=solve_system(alpha0,n_nodes,init_coord,approx_period);
        % numerically solves the problem
        
        
        % alpha1.value{2}=[-2*pi*1i,mu1,-mu1,-1]; 
        % xBarGen1=solve_system(alpha1,m,init_coord,approx_period);
        % same for x1
        
        coefs_linear{3}(1)=0;
        for i=1:alpha0.size_scalar
            coefs_linear{3}(i)=-coefs_linear{1}(i,:)*xBarGen0.scalar.'-...
                sum(sum(squeeze(coefs_linear{2}(i,:,:)).*xBarGen0.vector,1));
        end
    otherwise 
        
        s=sprintf('The system %s was not found, program interrupted.\n',Name_system);
        error(s)
end


% if requested, the function will plot the first computed solution
if flag_plot
    scrsz = get(0,'ScreenSize');
    figure('Position',[1 scrsz(4)/2 scrsz(3)*3/4 scrsz(4)/2])
    if xBarGen0.size_vector==3
        %plot3(xBarGen1,'v');
        %hold on
        plot3(xBarGen0,'r');
    elseif xBarGen0.size_vector==2
        %plot2(xBarGen1,'v');
        %hold on
        plot2(xBarGen0,'r');
    else
        %plot(xBarGen1,'v');
        %hold on
        plot(xBarGen0,'r');
    end
    hold on
end

% NEWTON on the smaller system (parameter fixed)
[xBar0]=Newton_Xi(xBarGen0,alpha0,coefs_linear,maxiter,min_res);
coefs_linearshort=coefs_linear;
% short solution stored (parameter excluded)
xShort=xBar0;

% parameter is added
xBar0.size_scalar=xBar0.size_scalar+1;
xBar0.scalar=[xBar0.scalar,new_var];

% non-squared deivative is computed
coefs_linear{1}(:,end+1)=0;
full_DH0=Function_derivative(xBar0,alpha_coef,coefs_linear,0);
% flag 0 to have non-square output

% nullspace of the derivative is computed
x0_dot=null(full_DH0);

% el=x0_dot( alpha_coef.size_scalar+ xBar0.nodes+1);
% x0_dot=x0_dot *conj(el);
% % symmetric = @(x) (x+conj(x(end:-1:1)))/2;
% % 
% % x0_dot(1:alpha_coef.size_scalar)= real(x0_dot(1:alpha_coef.size_scalar));
% % for i=1:alpha_coef.size_vector
% %     indeces= alpha_coef.size_scalar+ (i-1)*(2*xBar0.nodes+1) +...
% %     1:(2*xBar0.nodes+1);
% %     x0_dot(indeces) = symmetric(x0_dot(indeces));
% % end
% x0_dot=x0_dot/norm(x0_dot);



x0_dot= Xi_vec2vec(symmetrise(vec2Xi_vec(x0_dot,xBar0.size_scalar,xBar0.size_vector,...
    xBar0.nodes)));
x0_dot=x0_dot/norm(x0_dot);


% outputs renamed
x=xBar0;x_dot=x0_dot; DH=full_DH0;
return
end





% function xdot_sym=symmetrise_xdot(xdot,x)
% xdot_sym= xdot;
% 
% length_nodes=2*x.nodes+1;
% sym_vec=@(v) ( v + conj(v(end:-1:1,:)))/2;
% 
% for j=1:x.size_vector
%     vector=x.size_scalar+(j-1)*length_nodes+1:x.size_scalar+j*length_nodes;
%     xdot_sym(vector)=sym_vec(xdot(vector));
% end
% end
