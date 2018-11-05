% some examples on van der Pol to run the code and check it out
% as of 31 August 2017: it runs

global azabaza
global nu
global use_intlab
global talkative
global RAD_MAX
global Display
Display = 1;
talkative = 1;
use_intlab = 0;
RAD_MAX = 10^-2;

if isempty(azabaza)
    addpath(genpath('../'))
    startintlab
    azabaza =1;
end


%% standard Van der Pol
%
% nu = 1.1;
% mu = 1.3;
% n_nodes = 50;
%
%
% string_van_der_pol = '- dot x1 + l1 x2 \n - dot x2 + mu l1 x2 - mu l1 x1^2 x2 - l1 x1'; % general van der pol
% string_van_der_pol_mu = strrep(string_van_der_pol, 'mu', num2str(mu)); % plugging in mu
%
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
% scal_eq = default_scalar_eq(sol);
% polynomial = from_string_to_polynomial_coef(string_van_der_pol_mu);
% F = full_problem(scal_eq, polynomial);
%
% % NEWTON
% [sol2,yBar,res,DFm,RES] =Newton_2(sol,F,30,10^-7);
%
%
% % validation
% DF =  derivative(F,sol,0);
% DF_mat = derivative_to_matrix(DF);
% A  = inv(DF_mat);
%
%
% use_intlab = 1;
%
% Y_vector = Y_bound_new(A,sol2,F);
% Z0_vector=Z0_bound(DF_mat,A,sol2);
% Z1_vector=Z1_bound_new(A,sol2,F);    % problem here
% Z2_vector= Z2_bound_new(A,sol2,F);
% [Imin,Imax]=find_negative(Z2_vector,Z1_vector,Z0_vector,Y_vector);
%
% % check
% if Imax>RAD_MAX
%     Imax = RAD_MAX;
% end
%
% if talkative
%     fprintf('\n The interval is [ %e, %e ].\n\n',Imin,Imax);
% end



%% Van der Pol tests
% 
% mu_vector = 1.1;%[0.01,0.5:0.5:3];
% nu_vector = [1.001,1.02,1.05];%,1.07,1.1,1.13,1.15,1.2];
% n_nodes_vector = 20:5:350;
% figure (1);% figure(2); hold on; figure(3); hold on;
% axes('FontSize',15)
% title('Interval of validation depending on nodes','Interpreter','Latex','FontSize',20);
% ylabel('Interval','Interpreter','Latex','FontSize',20)
% xlabel('nodes','Interpreter','Latex','FontSize',20)
% % %set(gca, 'XTick', nu_vector, 'XTickLabel', {'1.001','1.02','1.05','1.07','1.1','1.13','1.15','1.2'});
% hold on
% set(gcf,'Visible', 'off'); 
% set(gca,'yscale','log');
% colormap winter
% colors = colormap;
% 
% string_van_der_pol = '- dot x1 + l1 x2 \n - dot x2 + mu l1 x2 - mu l1 x1^2 x2 - l1 x1'; % general van der pol
% n_nodes = 50;
% % construction of approximate solution (taking the circle every time, easy)
%         sin_four = zeros(1,2*n_nodes+1);
%         cos_four = sin_four;
%         cos_four(n_nodes) = 1/2; cos_four(n_nodes+2) = 1/2;
%         sin_four(n_nodes) = 1i/2; sin_four(n_nodes+2) = -1i/2;
%         x1 = sin_four;
%         x2 = cos_four;
%         
%         % constructing the problem, vector field and simple phase condition
%         sol2 = Xi_vector(1, [x1;x2]);% + 0.01*rand(2,1+2*n_nodes));
% 
% nu = nu_vector(1);
% for index_mu =1:length(mu_vector)
%     mu = mu_vector(index_mu);
%     string_van_der_pol_mu = strrep(string_van_der_pol, 'mu', num2str(mu)); % plugging in mu
%     polynomial = from_string_to_polynomial_coef(string_van_der_pol_mu);
%     
%     
%     %n_nodes = 50;
%     for n_nodes = n_nodes_vector
%         use_intlab = 0;
%         
%         sol2 = reshape_Xi(sol2,n_nodes);
%         scal_eq = default_scalar_eq(sol2);
%         F = full_problem(scal_eq, polynomial);
%         
%         % NEWTON
%         nu = nu_vector(1);
%         try
%             [sol2,yBar,res,DFm,RES] =Newton_2(sol2,F,30,10^-7);
%             %plot2(sol2, 'color',colors(8*index_mu,:),'LineWidth',2); hold on; %break
%         catch
%             break
%         end
%             
%             % validation
%             DF =  derivative(F,sol2,0);
%             DF_mat = derivative_to_matrix(DF);
%             A  = inv(DF_mat);
%             
%             
%             use_intlab = 1;
%             %nu = nu_vector;
%             for nu_index = 1:length(nu_vector)
%                 nu = nu_vector(nu_index);
%                 try
%                     Y_vector = Y_bound_new(A,sol2,F);
%                     Z0_vector=Z0_bound(DF_mat,A,sol2);
%                     Z1_vector=Z1_bound_new(A,sol2,F);
%                     Z2_vector= Z2_bound_new(A,sol2,F);
%                     [Imin,Imax]=find_negative(Z2_vector,Z1_vector,Z0_vector,Y_vector);
%                     
%                     % check
%                     if Imax>RAD_MAX
%                         Imax = RAD_MAX;
%                     end
%
%                     if talkative
%                         fprintf('\n The interval is [ %e, %e ].\n\n',Imin,Imax);
%                     end
%
%                     figure(1)
%                     semilogy(n_nodes, Imin, '.','color', colors(end - nu_index*8,:),'MarkerSize', 30);
%                     semilogy(n_nodes, Imax, '.','color', colors(nu_index*12,:),'MarkerSize', 30);
%                     % plot(nu, n_nodes, '.g','MarkerSize', 30)
%                     %k_min = n_nodes;
%                     %break
%                 catch
%                     %if n_nodes>100
%                     %    error('')
%                     %end
%                     %figure(index_mu)
%                     %plot(nu, n_nodes, '.r','MarkerSize', 30)
%                 end
%             end
%
%     end
% end
% legend('\nu=1.001, Imin','\nu=1.02, Imin','\nu=1.05, Imin','\nu=1.001, Imax','\nu=1.02, Imax','\nu=1.05, Imax','interpreter','Latex');
% figure(1)

%% weird van der Pol with powers
%%% Van der Pol tests


mu_vector = 1;%[0.01,0.5:0.5:3];
nu_vector = [1.001];%,1.02,1.05];%,1.07,1.1,1.13,1.15,1.2];
n_nodes_vector = 100;%20:5:350;
power_vector = 2:1:30;

set(0, 'DefaultFigureVisible', 'off');

figure (1); %hold on; figure(2); hold on; figure(3); hold on;
axes('FontSize',15)
title('Interval of validation depending on order of perturbation','Interpreter','Latex','FontSize',20);
ylabel('Interval','Interpreter','Latex','FontSize',20)
xlabel('power','Interpreter','Latex','FontSize',20)
set(gca,'yscale','log');
% % %set(gca, 'XTick', nu_vector, 'XTickLabel', {'1.001','1.02','1.05','1.07','1.1','1.13','1.15','1.2'});
hold on
set(gcf,'Visible', 'off');

figure (2);% figure(2); hold on; figure(3); hold on;
axes('FontSize',15)
title('Time of validation depending on order of perturbation','Interpreter','Latex','FontSize',20);
ylabel('time (s)','Interpreter','Latex','FontSize',20)
xlabel('power','Interpreter','Latex','FontSize',20)
hold on
set(gcf,'Visible', 'off');

figure(3);
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
        
        figure(3)
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
                
                figure(1)
                semilogy(power, Imin, '.','color', colors(end - nu_index*8,:),'MarkerSize', 30);
                semilogy(power, Imax, '.','color', colors(nu_index*12,:),'MarkerSize', 30);
                
                int_min(power)=Imin;
                int_max(power)= Imax;
                
                figure(2)
                plot(power, T,'.','color', colors(nu_index*12,:),'MarkerSize', 30)
                
                times(power) = T;
                % plot(nu, n_nodes, '.g','MarkerSize', 30)
                %k_min = n_nodes;
                %break
                %fprintf(2,[ sprintf('succeeded for power = %d',power) char(10)])
                %pause(2)
            catch
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

figure(1); 
saveas(gcf,'saved elements/interval_on_order','epsc');
saveas(gcf,'saved elements/interval_on_order','fig');
figure(2); 
saveas(gcf,'saved elements/time_on_order','epsc');
saveas(gcf,'saved elements/time_on_order','fig');
figure(3); 
saveas(gcf,'saved elements/solution_on_order','epsc');
saveas(gcf,'saved elements/solution_on_order','fig');
save('saved elements/results_power','int_min','int_max','times')
close all

return

%% weird van der Pol of higher dimension
%%% Van der Pol tests

mu_vector = 1;
nu_vector = [1.001];
n_nodes_vector = 50;
dimension_vector = 1:10;%[10:1:30];

set(0, 'DefaultFigureVisible', 'off');

figure (1); %figure(2); hold on; figure(3); hold on;
axes('FontSize',15)
title('Interval of validation depending on dimension of the system','Interpreter','Latex','FontSize',20);
ylabel('Interval','Interpreter','Latex','FontSize',20)
xlabel('dimension','Interpreter','Latex','FontSize',20)
set(gca,'yscale','log');
% % %set(gca, 'XTick', nu_vector, 'XTickLabel', {'1.001','1.02','1.05','1.07','1.1','1.13','1.15','1.2'});
hold on
% %set(gcf,'Visible', 'off');
% 
% figure (2);% figure(2); hold on; figure(3); hold on;
% axes('FontSize',15)
% title('Time of validation depending on dimension of the system','Interpreter','Latex','FontSize',20);
% ylabel('time (s)','Interpreter','Latex','FontSize',20)
% xlabel('dimension','Interpreter','Latex','FontSize',20)
% hold on
%set(gcf,'Visible', 'off');
% 
% figure(3);
% axes('FontSize',15)
% title('Solutions','Interpreter','Latex','FontSize',20);
% ylabel('y','Interpreter','Latex','FontSize',20)
% xlabel('x','Interpreter','Latex','FontSize',20)
% hold on 
% set(gcf,'Visible', 'off');

% set(gcf,'Visible', 'off');
% set(gca,'yscale','log');
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
        
        %figure(3)
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
                
                figure(1)
                semilogy(dimension, Imin, '.','color', colors(end - nu_index*8,:),'MarkerSize', 30);
                semilogy(dimension, Imax, '.','color', colors(nu_index*12,:),'MarkerSize', 30);
                
                int_min(dimension) = Imin;
                int_max(dimension) = Imax;
                
                figure(2)
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

% figure(1); 
% saveas(gcf,'saved elements/interval_on_dim','epsc');
% saveas(gcf,'saved elements/interval_on_dim','fig');
% figure(2); 
% saveas(gcf,'saved elements/time_on_dim','epsc');
% saveas(gcf,'saved elements/time_on_dim','fig');
% figure(3); 
% saveas(gcf,'saved elements/solution_on_dim','epsc');
% saveas(gcf,'saved elements/solution_on_dim','fig');
% save('saved elements/results_dim','int_min','int_max','times')

