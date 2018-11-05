% turned into a script
% function alpha=interactive_constructor(Name_system)
% function alpha=interactive_constructor()
%
% This function help constructing the element saving the coefficients of
% the ODE system under consideration. It will save the constructed
% coefficient to  file, in order to retreive them easily at the next
% iteration.
%

format compact

do_all=1;
if ~exist('Name_system','var')
    Name_system=input('Insert the name of the system to be solved: ','s');
end

try
    load(Name_system)
    if ~exist('alpha','var') && ~exist('alpha_coef','var')
        error('In the file needs to be a coefs variable named alpha');
    end
    if ~exist('alpha_coef','var')
        alpha_coef=alpha;
    end
    sflag=input('The system is already saved, do you want to overwrite it?','s');
    if sflag~='Y' || sflag~=''
        do_all=0;
    end
catch
    fprintf('The give system is unknown, we will continue with the procedure of constructing \nthe coefficients\n\n');
end

fprintf('\n')

if do_all
    flag=true;
    while flag
        size_scalar=input('Number of scalar unknowns: ');
        if int8(size_scalar)~=size_scalar
            disp('requires an integer input');
        else
            flag=false;
        end
    end
    if size_scalar>1
        warning('this feature has not yet been implemented')
    end
    
    flag=true;
    while flag
        size_vector=input('Number of function unknowns: ');
        if int8(size_vector)~=size_vector
            disp('requires an integer input');
        else
            flag=false;
        end
    end
    
    Koef=cell(size_vector,1);
    powers=cell(size_vector,1);
    powers_scalar=cell(size_vector,1);
    non_zero=cell(size_vector,1);
    deg_scal=0;
    deg_vec=0;
    fprintf('The number of non-zero coefficients for the vector functions\nwill be asked soon. Consider that the first of this coefficients \nis always connected to a derivative over time. \n')
    
    for i=1:size_vector
        flag=true;
        while flag
            fprintf('Considering the function %d\n',i);
            non_zero{i}=input('Insert the number of non-zero cefficients: ');
            if int8(non_zero{i})~=non_zero{i} && non_zero{i}<=0
                disp('requires a strictly positive integer input');
            else
                flag=false;
            end
        end
        Koef{i}=zeros(non_zero{i},1);
        powers{i}=zeros(size_vector,non_zero{i});
        powers_scalar{i}=zeros(size_scalar,non_zero{i});
        for j=1:non_zero{i}
            flag=true;
            while flag
                fprintf('For term %d, input the ',j);
                Koef{i}(j)=input('coefficient:');
                if ~isfloat(Koef{i}(j))
                    disp('requires a scalar input');
                else
                    flag=false;
                end
            end
            
            for k=1:size_scalar
                flag=true;
                while flag
                    fprintf('and the power of the %d-th scalar: ',k);
                    powers_scalar{i}(k,j)=input('');
                    if int8(powers_scalar{i}(k,j))~=powers_scalar{i}(k,j) && powers_scalar{i}(k,j)<=0
                        disp('requires a strictly positive integer input');
                    else
                        flag=false;
                    end
                end
            end
            deg_scal=max(deg_scal,sum(powers_scalar{i}(:,j)));
            for k=1:size_vector
                flag=true;
                while flag
                    fprintf('and the power of the %d-th function: ',k);
                    powers{i}(k,j)=input('');
                    if int8(powers{i}(k,j))~=powers{i}(k,j) && powers{i}(k,j)<=0
                        disp('requires a strictly positive integer input');
                    else
                        flag=false;
                    end
                end
            end
            deg_vec=max(deg_vec,sum(powers{i}(:,j)));
        end
    end
    
    alpha_coef=coefs(size_scalar,size_vector,deg_scal,deg_vec,non_zero,powers_scalar,powers,Koef);
    
    sflag=input('Do you want to save this system under the name given at the beginning?','s');
    if strcmp(sflag,'') || strcmp(sflag,'Y')|| strcmp(sflag,'y')
        save(Name_system,'alpha_coef');
    end
    
    
end

% % create also the coefficients for the linear part
% we have
%    g(x_1,x_2,...,x_M,^a_1,^a_2,...,^a_N) \in R^M
% more precisely
%    g()_l = sum_(i=1)^M g^1_{li} x_i + ...
%               sum_{i=1}^N\sum_{k=-m}^m g^2_{lik} ^a_{ik}
% with g^1 \in R^{M x M} and g^2 \in R^{M x N x (2m+1) } 
%
% here we create g^1 and g^2 calling them coefs_linear{1} and {2}
% respectively

coefs_linear=cell(3,1);
coefs_linear{1}=zeros(alpha_coef.size_scalar,1);
coefs_linear{2}=zeros(alpha_coef.size_scalar,alpha_coef.size_vector,2*n_nodes+1);
coefs_linear{3}=zeros(alpha_coef.size_scalar,1);% This is the vector of constants

% in the standard case (NO continuation yet), all the coefficients are 0
% exept the one refering to the coefficients from -k0 to k0 of the first
% sequence
k0=0;
coefs_linear{2}(1,1,1+k0:end-k0)=1; % MOST SIMPLE CASE

%alpha.linear=coefs_linear;



