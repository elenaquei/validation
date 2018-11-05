% some examples on van der Pol to run the code and check it out
% as of 31 August 2017: it runs but not hte last averaging part

%%
global nu
global use_intlab
global talkative
global RAD_MAX
global Display
Display = 1;
talkative = 1;
use_intlab = 0;
RAD_MAX = 10^-1;

try
    intval(1);
catch
    addpath(genpath('./'))
    startintlab
end
SAVE_IT_UP = 1;
%path_to_save = './saved elements/article1/';
path_to_save = './saved elements/dump/';

CASE=4;
%% Van der Pol tests - bounds depending on number of nodes and nu
% a pretty big mu is necessary, otherwise the solution is too simple
if CASE==1 %|| CASE ==0
    mu_vector = 1.1;
    nu_vector = 1:0.0125:1.2;% 1:0.0125:1.2;
    n_nodes_vector = 25:5:170;
    
    interval_min = zeros(length(nu_vector),length(n_nodes_vector));
    interval_max = zeros(length(nu_vector),length(n_nodes_vector));
    
    h_2 = figure;
    set(h_2,'Visible', 'on');
    axes('FontSize',15)
    title('Success or failure on nodes and $\nu$','Interpreter','Latex','FontSize',20);
    ylabel('Interval','Interpreter','Latex','FontSize',20)
    xlabel('nodes','Interpreter','Latex','FontSize',20)
    hold on;
    h_3=figure;
    set(h_3,'Visible', 'on');
    axes('FontSize',15)
    set(gca,'yscale','log');
    title('??','Interpreter','Latex','FontSize',20);
    ylabel('Interval','Interpreter','Latex','FontSize',20)
    xlabel('nodes','Interpreter','Latex','FontSize',20)
    hold on;
    % %set(gca, 'XTick', nu_vector, 'XTickLabel', {'1.001','1.02','1.05','1.07','1.1','1.13','1.15','1.2'});
    % hold on
    colormap hsv
    colors = colormap;
    
    colormap spring
    colors2 = colormap;
    
    string_van_der_pol = '- dot x1 + l1 x2 \n - dot x2 + mu l1 x2 - mu l1 x1^2 x2 - l1 x1'; % general van der pol
    n_nodes = 60;
    % construction of approximate solution (taking the circle every time, easy)
    sin_four = zeros(1,2*n_nodes+1);
    cos_four = sin_four;
    cos_four(n_nodes) = 1/2; cos_four(n_nodes+2) = 1/2;
    sin_four(n_nodes) = 1i/2; sin_four(n_nodes+2) = -1i/2;
    x1 = sin_four;
    x2 = cos_four;
    
    % constructing the problem, vector field and simple phase condition
    sol2 = Xi_vector(1, [x1;x2]);% + 0.01*rand(2,1+2*n_nodes));
    nu = nu_vector(1);
    
    %nu = nu_vector(1);
    for index_mu =1:1%length(mu_vector)
        mu = mu_vector(index_mu);
        string_van_der_pol_mu = strrep(string_van_der_pol, 'mu', num2str(mu)); % plugging in mu
        polynomial = from_string_to_polynomial_coef(string_van_der_pol_mu);
        scal_eq = default_scalar_eq(sol2);
        F = full_problem(scal_eq, polynomial);
        [sol2] =Newton_2(sol2,F,30,10^-7);
        
        %n_nodes = 50;
        for i_n_nodes = 1:length(n_nodes_vector)
            n_nodes = n_nodes_vector(i_n_nodes);
            use_intlab = 0;
            
            %if n_nodes~=sol2_test.nodes
            sol2 = reshape_Xi(sol2,n_nodes);
            scal_eq = default_scalar_eq(sol2);
            F = full_problem(scal_eq, polynomial);
            
            % NEWTON
            nu = nu_vector(1);
            try
                [sol2,yBar,res,DFm,RES] =Newton_2(sol2,F,30,10^-7);
                %plot2(sol2, 'color',colors(8*index_mu,:),'LineWidth',2); hold on; %break
            catch
                continue
            end
            %else
            %    sol2 = sol2_test;
            %    sol2_newton = sol2;
            %end
            
            % validation
            DF =  derivative(F,sol2,0);
            DF_mat = derivative_to_matrix(DF);
            cond_normal(i_n_nodes) = cond(DF_mat);
            A  = inv(DF_mat);
            
            DF_newton =  derivative(F,sol2,0);
            DF_mat_newton = derivative_to_matrix(DF_newton);
            cond_newton(i_n_nodes) = cond(DF_mat_newton);
            A_newton  = inv(DF_mat_newton);
            
            use_intlab = 1;
            %nu = nu_vector;
            for nu_index = 1:length(nu_vector)
                nu = nu_vector(nu_index);
                try
                    %                 Y_vector = Y_bound_new(A,sol2,F);
                    %                 Z0_vector=Z0_bound(DF_mat,A,sol2);
                    %                 Z1_vector=Z1_bound_new(A,sol2,F);
                    %                 Z2_vector= Z2_bound_new(A,sol2,F);
                    %                 [Imin,Imax]=find_negative(Z2_vector,Z1_vector,Z0_vector,Y_vector);
                    
                    
                    Y_vector = Y_bound_new(A_newton,sol2,F);
                    Z0_vector=Z0_bound(DF_mat_newton,A_newton,sol2);
                    Z1_vector=Z1_bound_new(A_newton,sol2,F);
                    Z2_vector= Z2_bound_new(A_newton,sol2,F);
                    [Imin_newton,Imax_newton]=find_negative(Z2_vector,Z1_vector,Z0_vector,Y_vector);
                    
                    % check
                    %if Imax>RAD_MAX
                    %    Imax = RAD_MAX;
                    %end
                    
                    
                    if Imax_newton>RAD_MAX
                        Imax_newton = RAD_MAX;
                    end
                    
                    if talkative
                        fprintf('\n The interval is [ %e, %e ].\n\n',Imin_newton,Imax_newton);
                    end
                    interval_min(nu_index,i_n_nodes) = Imin_newton;
                    interval_max(nu_index,i_n_nodes) = Imax_newton;
                    %figure(2)
                    %semilogy(n_nodes, Imin, '.','color', colors(end - nu_index*8,:),'MarkerSize', 30);
                    %semilogy(n_nodes, Imax, '.','color', colors(nu_index*12,:),'MarkerSize', 30);
                    %                 if n_nodes ==45
                    %                     x_imax(nu_index)=semilogy(n_nodes, Imin_newton, '.','color', colors(end - nu_index*8,:),'MarkerSize', 30);
                    %                     x_imin(nu_index)=semilogy(n_nodes, Imax_newton, '.','color', colors(nu_index*4,:),'MarkerSize', 30);
                    %                 else
                    %if mod(nu_index,4)==0
                    figure(h_3)
                    semilogy(n_nodes, Imin_newton, '.','color', colors(end - nu_index*2,:),'MarkerSize', 30);
                    semilogy(n_nodes, Imax_newton, '.','color', colors(nu_index,:),'MarkerSize', 30);
                    %end
                    figure(h_2)
                    plot(nu, n_nodes, '.g','MarkerSize', 30)
                    %k_min = n_nodes;
                    %break
                catch e
                    %if n_nodes>100
                    %    error('')
                    %end
                    figure(h_2)
                    plot(nu, n_nodes, '.r','MarkerSize', 30)
                end
            end
            
        end
    end
    %legend([x_imin,x_imax],'\nu=1.001, Imin','\nu=1.02, Imin','\nu=1.05, Imin','\nu=1.001, Imax','\nu=1.02, Imax','\nu=1.05, Imax','location','westoutside')
    % legend('\nu=1.001, Imin','\nu=1.02, Imin','\nu=1.05, Imin','\nu=1.001, Imax','\nu=1.02, Imax','\nu=1.05, Imax','interpreter','Latex');
    if SAVE_IT_UP
        figure(h_2)
        set(gcf, 'Position', get(0, 'Screensize'));
        saveas(gcf,strcat(path_to_save,'success_on_mu_and_nu'),'epsc');
        saveas(gcf,strcat(path_to_save,'success_on_mu_and_nu'),'fig');
        figure(h_3)
        set(gcf, 'Position', get(0, 'Screensize'));
        saveas(gcf,strcat(path_to_save,'interval_on_mu_and_nu'),'epsc');
        saveas(gcf,strcat(path_to_save,'interval_on_mu_and_nu'),'fig');
        save(strcat(path_to_save,'interval_on_mu_and_nu'), 'interval_min','interval_max');
    end
    
end

%% Van der Pol tests - minimum number of nodes depending on mu
% starting with mu really small, we go up adding nodes if necessary to get
% the validation running
if CASE==2 %|| CASE ==0
    mu_vector = 0.4:0.2:7;
    nu = 1.001;
    min_nodes = 10;
    %n_nodes_vector = 10;%:5:500;
    h_4 = figure;
    axes('FontSize',15)
    title('Minimum number of nodes depending on $\mu$ ','Interpreter','Latex','FontSize',20);
    ylabel('nodes','Interpreter','Latex','FontSize',20)
    xlabel('$\mu$','Interpreter','Latex','FontSize',20)
    hold on
    set(h_4,'Visible', 'off');
    colormap winter
    colors = colormap;
    
    string_van_der_pol = '- dot x1 + l1 x2 \n - dot x2 + mu l1 x2 - mu l1 x1^2 x2 - l1 x1'; % general van der pol
    n_nodes = min_nodes;
    % construction of approximate solution (taking the circle every time, easy)
    sin_four = zeros(1,2*n_nodes+1);
    cos_four = sin_four;
    cos_four(n_nodes) = 1/2; cos_four(n_nodes+2) = 1/2;
    sin_four(n_nodes) = 1i/2; sin_four(n_nodes+2) = -1i/2;
    x1 = sin_four;
    x2 = cos_four;
    
    % constructing the problem, vector field and simple phase condition
    sol2 = Xi_vector(1, [x1;x2]);% + 0.01*rand(2,1+2*n_nodes));
    %nu = nu_vector(1);
    
    for index_mu =1:length(mu_vector)
        use_intlab=0;
        mu = mu_vector(index_mu);
        string_van_der_pol_mu = strrep(string_van_der_pol, 'mu', num2str(mu)); % plugging in mu
        polynomial = from_string_to_polynomial_coef(string_van_der_pol_mu);
        scal_eq = default_scalar_eq(sol2);
        F = full_problem(scal_eq, polynomial);
        [sol2] =Newton_2(sol2,F,30,10^-7);
        n_nodes = sol2.nodes;
        validation = 0;
        tries = 0;
        while ~validation && tries <40
            tries = tries+1;
            use_intlab = 0;
            
            
            sol2 = reshape_Xi(sol2,n_nodes);
            scal_eq = default_scalar_eq(sol2);
            F = full_problem(scal_eq, polynomial);
            
            % NEWTON
            nu = nu_vector(1);
            try
                [sol2,yBar,res,DFm,RES] =Newton_2(sol2,F,30,10^-7);
                %plot2(sol2, 'color',colors(8*index_mu,:),'LineWidth',2); hold on; %break
            catch
                n_nodes = n_nodes +1;
                continue
            end
            
            
            % validation
            DF =  derivative(F,sol2,0);
            DF_mat = derivative_to_matrix(DF);
            A  = inv(DF_mat);
            
            use_intlab = 1;
            try
                Y_vector = Y_bound_new(A,sol2,F);
                Z0_vector=Z0_bound(DF_mat,A,sol2);
                Z1_vector=Z1_bound_new(A,sol2,F);
                Z2_vector= Z2_bound_new(A,sol2,F);
                [Imin,Imax]=find_negative(Z2_vector,Z1_vector,Z0_vector,Y_vector);
                
                
                % check
                if Imax>RAD_MAX
                    Imax = RAD_MAX;
                end
                
                if talkative
                    fprintf('\n The interval is [ %e, %e ].\n\n',Imin,Imax);
                end
                
                figure(h_4)
                plot(mu, n_nodes, '.','color', colors(5,:),'MarkerSize', 30)
                mu_on_nodes(index_mu,1:2) = [mu,n_nodes];
                %k_min = n_nodes;
                %break
                validation = 1;
            catch
                n_nodes = n_nodes +1;
            end
            use_intlab =0;
        end
        
    end
    if SAVE_IT_UP
        figure(h_4)
        set(gcf, 'Position', get(0, 'Screensize'));
        saveas(gcf,strcat(path_to_save,'min_nodes_on_mu'),'epsc');
        saveas(gcf,strcat(path_to_save,'min_nodes_on_mu'),'fig');
        save(strcat(path_to_save,'min_nodes_on_mu'), 'mu_on_nodes');
    end
    %legend('\nu=1.001, Imin','\nu=1.02, Imin','\nu=1.05, Imin','\nu=1.001, Imax','\nu=1.02, Imax','\nu=1.05, Imax','location','westoutside')
    %legend('\nu=1.001, Imin','\nu=1.02, Imin','\nu=1.05, Imin','\nu=1.001, Imax','\nu=1.02, Imax','\nu=1.05, Imax','interpreter','Latex');
    %figure(3)
    
end
%% weird van der Pol with powers
%%% Van der Pol tests

%close all

if CASE ==3 || CASE ==0
    mu_vector = 1;%[0.01,0.5:0.5:3];
    nu_vector = 1.001;%,1.02,1.05];%,1.07,1.1,1.13,1.15,1.2];
    n_nodes_vector = 100;%20:5:350;
    power_vector = 2:1:30;
    
    
    
    h_5=figure; %hold on; figure(2); hold on; figure(3); hold on;
    axes('FontSize',15)
    title('Interval of validation depending on order of perturbation','Interpreter','Latex','FontSize',20);
    ylabel('Interval','Interpreter','Latex','FontSize',20)
    xlabel('power','Interpreter','Latex','FontSize',20)
    set(gca,'yscale','log');
    % % %set(gca, 'XTick', nu_vector, 'XTickLabel', {'1.001','1.02','1.05','1.07','1.1','1.13','1.15','1.2'});
    hold on
    set(gcf,'Visible', 'off');
    
    h_6=figure;% figure(2); hold on; figure(3); hold on;
    axes('FontSize',15)
    title('Time of validation depending on order of perturbation','Interpreter','Latex','FontSize',20);
    ylabel('time (s)','Interpreter','Latex','FontSize',20)
    xlabel('power','Interpreter','Latex','FontSize',20)
    hold on
    set(gcf,'Visible', 'off');
    
    h_7=figure;
    axes('FontSize',15)
    title('Solutions','Interpreter','Latex','FontSize',20);
    ylabel('y','Interpreter','Latex','FontSize',20)
    xlabel('x','Interpreter','Latex','FontSize',20)
    hold on
    set(gcf,'Visible', 'off');
    
    % set(gcf,'Visible', 'off');
    % set(gca,'yscale','log');
    colormap winter
    colors = colormap;
    
    times = zeros(length(power_vector),1);
    int_min = times;
    int_max = times;
    Z2_vector=zeros(length(power_vector),3);
    Z1_vector=zeros(length(power_vector),3);
    Z0_vector=zeros(length(power_vector),3);
    Y_vector=zeros(length(power_vector),3);
    
    
    string_van_der_pol = '- dot x1 + l1 x2 \n - dot x2 + mu l1 x2 - mu l1 x1^2 x2 - l1 x1+eps x1 ^pow'; % general van der pol
    n_nodes = 50;
    % construction of approximate solution (taking the circle every time, easy)
    sin_four = zeros(1,2*n_nodes+1);
    cos_four = sin_four;
    cos_four(n_nodes) = 1/2; cos_four(n_nodes+2) = 1/2;
    sin_four(n_nodes) = 1i/2; sin_four(n_nodes+2) = -1i/2;
    x1 = sin_four;
    x2 = cos_four;
    
    % constructing the problem, vector field and simple phase condition
    sol = Xi_vector(1, [x1;x2]);% + 0.01*rand(2,1+2*n_nodes));
    sol2 = sol;
    
    nu = nu_vector(1);
    mu = mu_vector(1);
    
    %times = zeros(max(power_vector));
    %int_min = times;
    %int_max = times;
    
    
    for index_pow =1:length(power_vector)
        power = power_vector(index_pow); % change this to change the power
        if strcmp(num2str(0.1/(2^power),'%1.9f'),num2str(0,'%1.9f'))
            break
        end
        string_van_der_pol_mu = strrep(string_van_der_pol, 'mu', num2str(mu)); % plugging in mu
        string_van_der_pol_mu = strrep(string_van_der_pol_mu, 'pow', num2str(power));
        string_van_der_pol_mu = strrep(string_van_der_pol_mu, 'eps', num2str(0.1/(2^power),'%1.9f'));
        polynomial = from_string_to_polynomial_coef(string_van_der_pol_mu);
        
        %n_nodes = 50;
        for n_nodes = n_nodes_vector
            use_intlab = 0;
            
            sol2 = reshape_Xi(sol2,n_nodes);
            scal_eq = default_scalar_eq(sol2);
            F = full_problem(scal_eq, polynomial);
            
            % NEWTON
            nu = nu_vector(1);
            try
                [sol2,yBar,res,DFm,RES] =Newton_2(sol2,F,30,10^-7);
                
            catch
                fprintf('Newton failed for pow = %i\n',power)
                break
            end
            
            figure(h_7)
            plot2(sol2, 'color',colors(2*power,:),'LineWidth',2);
            
            
            % validation
            DF =  derivative(F,sol2,0);
            DF_mat = derivative_to_matrix(DF);
            A  = inv(DF_mat);
            
            
            use_intlab = 1;
            
            for nu_index = 1:length(nu_vector)
                nu = nu_vector(nu_index);
                try
                    T1 = cputime;
                    Y_vector (index_pow,:)= Y_bound_new(A,sol2,F);
                    Z0_vector(index_pow,:)=Z0_bound(DF_mat,A,sol2);
                    Z1_vector(index_pow,:)=Z1_bound_new(A,sol2,F);
                    Z2_vector(index_pow,:)= Z2_bound_new(A,sol2,F);
                    [Imin,Imax]=find_negative(Z2_vector (index_pow,:),Z1_vector (index_pow,:),Z0_vector (index_pow,:),Y_vector (index_pow,:));
                    
                    % check
                    if Imax>RAD_MAX
                        Imax = RAD_MAX;
                    end
                    
                    if talkative
                        fprintf('\n The interval is [ %e, %e ].\n\n',Imin,Imax);
                    end
                    T = cputime - T1;
                    
                    figure(h_5)
                    semilogy(power, Imin, '.','color', colors(end - nu_index*8,:),'MarkerSize', 30);
                    semilogy(power, Imax, '.','color', colors(nu_index*12,:),'MarkerSize', 30);
                    
                    int_min(index_pow)=Imin;
                    int_max(index_pow)= Imax;
                    
                    figure(h_6)
                    plot(power, T,'.','color', colors(nu_index*12,:),'MarkerSize', 30)
                    
                    times(index_pow) = T;
                    % plot(nu, n_nodes, '.g','MarkerSize', 30)
                    %k_min = n_nodes;
                    %break
                    %fprintf(2,[ sprintf('succeeded for power = %d',power) char(10)])
                    %pause(2)
                catch e
                    fprintf(2,[ sprintf('failed for power = %d',power) char(10)])
                    % pause(2)
                    %if n_nodes>100
                    %    error('')
                    %end
                    %figure(index_mu)
                    %plot(nu, n_nodes, '.r','MarkerSize', 30)
                end
            end
            
        end
    end
    
    % figure(10); plot(norm_A)
    % figure(10)
    % figure(11); plot(omega)
    % figure(11)
    % title('Omega')
    % figure(10)
    % title('norm A')
    if SAVE_IT_UP
        figure(h_5)
        set(gcf, 'Position', get(0, 'Screensize'));
        saveas(gcf,strcat(path_to_save,'interval_on_order'),'epsc');
        saveas(gcf,strcat(path_to_save,'interval_on_order'),'fig');
        figure(h_6)
        
        set(gcf, 'Position', get(0, 'Screensize'));
        saveas(gcf,strcat(path_to_save,'time_on_order'),'epsc');
        saveas(gcf,strcat(path_to_save,'time_on_order'),'fig');
        
        figure(h_7);
        set(gcf, 'Position', get(0, 'Screensize'));
        saveas(gcf,strcat(path_to_save,'solution_on_order'),'epsc');
        saveas(gcf,strcat(path_to_save,'solution_on_order'),'fig');
        
        save(strcat(path_to_save,'results_power'),'int_min','int_max','times','power_vector')
    end
    %close all
    
   
end
%% weird van der Pol of higher dimension
%%% Van der Pol tests
if CASE == 4 || CASE ==0
    mu_vector = 1;
    nu_vector = 1.001;
    n_nodes_vector = 50;
    dimension_vector = 3:1:60;
    
    set(0, 'DefaultFigureVisible', 'off');
    
    h_8=figure; %figure(2); hold on; figure(3); hold on;
    axes('FontSize',15)
    title('Interval of validation depending on dimension of the system','Interpreter','Latex','FontSize',20);
    ylabel('Interval','Interpreter','Latex','FontSize',20)
    xlabel('dimension','Interpreter','Latex','FontSize',20)
    set(gca,'yscale','log');
    % % %set(gca, 'XTick', nu_vector, 'XTickLabel', {'1.001','1.02','1.05','1.07','1.1','1.13','1.15','1.2'});
    hold on
    % %set(gcf,'Visible', 'off');
    %
    h_9=figure;% figure(2); hold on; figure(3); hold on;
    axes('FontSize',15)
    title('Time of validation depending on dimension of the system','Interpreter','Latex','FontSize',20);
    ylabel('time (s)','Interpreter','Latex','FontSize',20)
    xlabel('dimension','Interpreter','Latex','FontSize',20)
    hold on
    set(gcf,'Visible', 'off');
    
     h_10=figure;
    axes('FontSize',15)
    title('Solutions','Interpreter','Latex','FontSize',20);
    ylabel('y','Interpreter','Latex','FontSize',20)
    xlabel('x','Interpreter','Latex','FontSize',20)
    hold on
    set(gcf,'Visible', 'off');
    
    set(gcf,'Visible', 'off');
    set(gca,'yscale','log');
    colormap winter
    colors = colormap;
    
    string_van_der_pol1 = '- dot x1 + l1 x2 \n';
    string_van_der_pol2 = '- dot xi + mu l1 xi - mu l1 x1^2 xi - l1 x1 \n'; % general van der pol
    n_nodes = 50;
    % construction of approximate solution (taking the circle every time, easy)
    sin_four = zeros(1,2*n_nodes+1);
    cos_four = sin_four;
    cos_four(n_nodes) = 1/2; cos_four(n_nodes+2) = 1/2;
    sin_four(n_nodes) = 1i/2; sin_four(n_nodes+2) = -1i/2;
    x1 = sin_four;
    x2 = cos_four;
    
    % constructing the problem, vector field and simple phase condition
    sol = Xi_vector(1, [x1;x2]);% + 0.01*rand(2,1+2*n_nodes));
    sol2 = sol;
    
    nu = nu_vector(1);
    mu = mu_vector(1);
    int_min = zeros(length(dimension_vector),1);
    int_max = zeros(length(dimension_vector),1);
    times = zeros(length(dimension_vector),1);
    for index_dim =1:length(dimension_vector)
        dimension = dimension_vector(index_dim); % change this to change the power
        string_2 = '';
        for i = 2:dimension
            string_temp = strrep(string_van_der_pol2, 'i',num2str(i));
            string_2=strcat(string_2,string_temp);
        end
        string_van_der_pol = strcat(string_van_der_pol1, string_2);
        string_van_der_pol_mu = strrep(string_van_der_pol, 'mu', num2str(mu)); % plugging in mu
        
        polynomial = from_string_to_polynomial_coef(string_van_der_pol_mu);
        
        %n_nodes = 50;
        for n_nodes = n_nodes_vector
            use_intlab = 0;
            
            sol2 = Xi_vector(1, [x1;repmat(x2,dimension-1,1)]);
            
            sol2 = reshape_Xi(sol2,n_nodes);
            scal_eq = default_scalar_eq(sol2);
            F = full_problem(scal_eq, polynomial);
            
            % NEWTON
            nu = nu_vector(1);
            try
                [sol2,yBar,res,DFm,RES] =Newton_2(sol2,F,30,10^-7);
                
            catch
                fprintf('Newton failed for dimension = %i\n',dimension)
                break
            end
            
            %figure(h_10)
            %plot2(sol2, 'color',colors(3*power,:),'LineWidth',2);
            
            
            % validation
            DF =  derivative(F,sol2,0);
            DF_mat = derivative_to_matrix(DF);
            A  = inv(DF_mat);
            
            
            use_intlab = 1;
            
            for nu_index = 1:length(nu_vector)
                nu = nu_vector(nu_index);
                try
                    T1 = cputime;
                    Y_vector = Y_bound_new(A,sol2,F);
                    Z0_vector=Z0_bound(DF_mat,A,sol2);
                    Z1_vector=Z1_bound_new(A,sol2,F);
                    Z2_vector= Z2_bound_new(A,sol2,F);
                    [Imin,Imax]=find_negative(Z2_vector,Z1_vector,Z0_vector,Y_vector);
                    
                    % check
                    if Imax>RAD_MAX
                        Imax = RAD_MAX;
                    end
                    
                    if talkative
                        fprintf('\n The interval is [ %e, %e ].\n\n',Imin,Imax);
                    end
                    T = cputime - T1;
                    
                    figure(h_8)
                    semilogy(dimension, Imin, '.','color', colors(end - nu_index*8,:),'MarkerSize', 30);
                    semilogy(dimension, Imax, '.','color', colors(nu_index*12,:),'MarkerSize', 30);
                    
                    int_min(dimension) = Imin;
                    int_max(dimension) = Imax;
                    
                    figure(h_9)
                    plot(dimension, T,'.','color', colors(nu_index*12,:),'MarkerSize', 30)
                    
                    times(dimension) = T;
                    % plot(nu, n_nodes, '.g','MarkerSize', 30)
                    %k_min = n_nodes;
                    %break
                    %fprintf(2,[ sprintf('succeeded for power = %d',power) char(10)])
                    %pause(2)
                catch
                    fprintf(2,[ sprintf('failed for dimension = %d',dimension) char(10)])
                    % pause(2)
                    %if n_nodes>100
                    %    error('')
                    %end
                    %figure(index_mu)
                    %plot(nu, n_nodes, '.r','MarkerSize', 30)
                end
            end
            
        end
    end
    if SAVE_IT_UP
        figure(h_8);
        set(gcf, 'Position', get(0, 'Screensize'));
        saveas(gcf,strcat(path_to_save,'interval_on_dim'),'epsc');
        saveas(gcf,strcat(path_to_save,'interval_on_dim'),'fig');
        figure(h_9);
        set(gcf, 'Position', get(0, 'Screensize'));
        saveas(gcf,strcat(path_to_save,'time_on_dim'),'epsc');
        saveas(gcf,strcat(path_to_save,'time_on_dim'),'fig');
        %figure(3);
        %set(gcf, 'Position', get(0, 'Screensize'));
        %saveas(gcf,'saved elements/solution_on_dim','epsc');
        %saveas(gcf,'saved elements/solution_on_dim','fig');
        save(strcat(path_to_save,'results_dim'),'int_min','int_max','times')
    end
    return
end


%% weird van der Pol of higher dimension with averaging
%%% Van der Pol tests  ILL CONDITIONED - nothing done
%
% mu_vector = 1;
% nu_vector = [1.001];
% n_nodes_vector = 50;
% dimension_vector = 1:60;
%
% %set(0, 'DefaultFigureVisible', 'off');
%
% figure (1); %figure(2); hold on; figure(3); hold on;
% axes('FontSize',15)
% title('Interval of validation depending on dimension of the system','Interpreter','Latex','FontSize',20);
% ylabel('Interval','Interpreter','Latex','FontSize',20)
% xlabel('dimension','Interpreter','Latex','FontSize',20)
% set(gca,'yscale','log');
% % % %set(gca, 'XTick', nu_vector, 'XTickLabel', {'1.001','1.02','1.05','1.07','1.1','1.13','1.15','1.2'});
% hold on
% % %set(gcf,'Visible', 'off');
% %
% figure (2);% figure(2); hold on; figure(3); hold on;
% axes('FontSize',15)
% title('Time of validation depending on dimension of the system','Interpreter','Latex','FontSize',20);
% ylabel('time (s)','Interpreter','Latex','FontSize',20)
% xlabel('dimension','Interpreter','Latex','FontSize',20)
% hold on
% set(gcf,'Visible', 'off');
%
% figure(3);
% axes('FontSize',15)
% title('Solutions','Interpreter','Latex','FontSize',20);
% ylabel('y','Interpreter','Latex','FontSize',20)
% xlabel('x','Interpreter','Latex','FontSize',20)
% hold on
% set(gcf,'Visible', 'off');
%
% set(gcf,'Visible', 'off');
% set(gca,'yscale','log');
% colormap winter
% colors = colormap;
%
% string_van_der_pol1 = '- dot xi + l1 Yi \n'; % Yi = y(i+1)
% string_van_der_pol2 = '- dot yi + mu l1 yi - l1 x1 +averaging \n'; % general van der pol
%
% string_xi = @(i) sprintf('x%i',2*i-1);% returns the string for xi
% string_yi = @(i) sprintf('x%i', 2*i); % returns the string for yi
% string_Yi = @(i,d) string_yi(mod(i,d)+1); % returns the string for Yi
% string_j_k = @(j,k, mu, D) sprintf('- %1.19f x%i x%i ', mu/D^2, 2*j-1,2*k-1); % return the string mu/D^2 * xj * xk
% string_j_j = @(j, mu, D) sprintf('- %1.19f x%i^2 ', mu/D^2, 2*j-1); % return the string mu/D^2 * xj ^2
%
% n_nodes = 50;
% % construction of approximate solution (taking the circle every time, easy)
% sin_four = zeros(1,2*n_nodes+1);
% cos_four = sin_four;
% cos_four(n_nodes) = 1/2; cos_four(n_nodes+2) = 1/2;
% sin_four(n_nodes) = 1i/2; sin_four(n_nodes+2) = -1i/2;
% x1 = sin_four;
% x2 = cos_four;
%
% % constructing the problem, vector field and simple phase condition
% sol = Xi_vector(1, [x1;x2]);% + 0.01*rand(2,1+2*n_nodes));
% sol2 = sol;
%
% nu = nu_vector(1);
% mu = mu_vector(1);
%
% for index_dim =1:length(dimension_vector)
%     dimension = dimension_vector(index_dim); % change this to change the power
%     averaging = '';
%     for j=1:dimension
%         for k =1:dimension
%             if k ~= j
%                 averaging = strcat(averaging, string_j_k(j,k,mu,dimension));
%             else
%                 averaging = strcat(averaging, string_j_j(j,mu,dimension));
%             end
%             yi = string_yi(i);
%             averaging = strcat(averaging, 'yi');
%         end
%     end
%
%     string_2 = '';
%     for i = 1:dimension
%         xi = string_xi(i);
%         yi = string_yi(i);
%         Yi = string_Yi(i,dimension);
%         string_x = strrep(strrep(string_van_der_pol1, 'xi',xi),'Yi',Yi);
%         string_y = strrep(string_van_der_pol2, '+averaging', averaging);
%         string_y = strrep(strrep(string_y, 'xi',xi),'yi',yi);
%         string_2=strcat(strcat(string_2,string_x),string_y);
%     end
%     %string_van_der_pol = strcat(string_van_der_pol1, string_2);
%     string_van_der_pol_mu = strrep(string_2, 'mu', num2str(mu)); % plugging in mu
%
%     polynomial = from_string_to_polynomial_coef(string_van_der_pol_mu);
%
%     %continue
%
%     %n_nodes = 50;
%     for n_nodes = n_nodes_vector
%         use_intlab = 0;
%
%         sol2 = Xi_vector(1, repmat([x1;x2],dimension,1));
%
%         sol2 = reshape_Xi(sol2,n_nodes);
%         scal_eq = default_scalar_eq(sol2);
%         F = full_problem(scal_eq, polynomial);
%
%         % NEWTON
%         nu = nu_vector(1);
%         try
%             [sol2,iter, yBar,res,DFm,RES] =Newton_2(sol2,F,30,10^-10);
%
%         catch
%             fprintf('Newton failed for dimension = %i\n',dimension)
%             break
%         end
%         if dimension ==1
%             sol_test =sol2;
%         end
%         %figure(3)
%         %plot2(sol2, 'color',colors(3*power,:),'LineWidth',2);
%
%
%         % validation
%         DF =  derivative(F,sol2,0);
%         DF_mat = derivative_to_matrix(DF);
%         A  = inv(DF_mat);
%
%
%         use_intlab = 1;
%
%         for nu_index = 1:length(nu_vector)
%             nu = nu_vector(nu_index);
%             try
%                 T1 = cputime;
%                 Y_vector = Y_bound_new(A,sol2,F);
%                 Z0_vector=Z0_bound(DF_mat,A,sol2);
%                 Z1_vector=Z1_bound_new(A,sol2,F);
%                 Z2_vector= Z2_bound_new(A,sol2,F);
%                 [Imin,Imax]=find_negative(Z2_vector,Z1_vector,Z0_vector,Y_vector);
%
%                 % check
%                 if Imax>RAD_MAX
%                     Imax = RAD_MAX;
%                 end
%
%                 if talkative
%                     fprintf('\n The interval is [ %e, %e ].\n\n',Imin,Imax);
%                 end
%                 T = cputime - T1;
%
%                 figure(1)
%                 semilogy(dimension, Imin, '.','color', colors(end - nu_index*8,:),'MarkerSize', 30);
%                 semilogy(dimension, Imax, '.','color', colors(nu_index*12,:),'MarkerSize', 30);
%
%                 int_min(dimension) = Imin;
%                 int_max(dimension) = Imax;
%
%                 figure(2)
%                 plot(dimension, T,'.','color', colors(nu_index*12,:),'MarkerSize', 30)
%
%                 times(dimension) = T;
%                 % plot(nu, n_nodes, '.g','MarkerSize', 30)
%                 %k_min = n_nodes;
%                 %break
%                 %fprintf(2,[ sprintf('succeeded for power = %d',power) char(10)])
%                 %pause(2)
%             catch
%                 fprintf(2,[ sprintf('failed for dimension = %d',dimension) char(10)])
%                 % pause(2)
%                 %if n_nodes>100
%                 %    error('')
%                 %end
%                 %figure(index_mu)
%                 %plot(nu, n_nodes, '.r','MarkerSize', 30)
%             end
%             use_intlab = 0;
%         end
%
%     end
% end
%
% figure(1);
% set(gcf, 'Position', get(0, 'Screensize'));
% saveas(gcf,'saved elements/interval_on_dim_averaging','epsc');
% saveas(gcf,'saved elements/interval_on_dim_averaging','fig');
% figure(2);
% set(gcf, 'Position', get(0, 'Screensize'));
% saveas(gcf,'saved elements/time_on_dim_averaging','epsc');
% saveas(gcf,'saved elements/time_on_dim_averaging','fig');
% figure(3);
% set(gcf, 'Position', get(0, 'Screensize'));
% saveas(gcf,'saved elements/solution_on_dim_averaging','epsc');
% saveas(gcf,'saved elements/solution_on_dim_averaging','fig');
% save('saved elements/results_dim_averaging','int_min','int_max','times')
%
