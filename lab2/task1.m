% ----- TASK 1 -----

clc, clear all, close all

disp('----- TASK 1 -----')
disp(' ')

color = {'b','r','g','k','m'};

% We're just to test whether we can sample from one of the densities. Thus
% we begin by defining two parameters. 

d_vek = 2:5; 

% This parameter represents the number of intervals we have, i.e. if we 
% have d = 3 we have two breakpoints.

beta = 1; % GOOD BETA IS 1

% Hyper parameter

rho_vek = [0.05 0.14 0.25 0.36];

% rho = 0.5 for d = 2
% rho = 0.75 for d = 3
% rho = 1 for d = 4
% rho = 1.25 for d = 5

% ---- THEORY -----
% http://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm
% -----------------


% From looking at the plots for different beta and rho, we can easily see
% that the MCMC algorithm is extremely robust for beta, wheras we can't
% vary rho that much.

% Tuning parameter

burn_in = 1e4;

% We define the burn in to be 5000

t_1 = 1851;t_end = 1963;N = 7.5e4 + burn_in;

% We then load the coal mine disasters

load coal_mine.mat

disa = coal_mine';

% We then plot the disasters 

[disaster_plot,years_plot] = hist(disa,20);
bar(years_plot,disaster_plot)
title('Histogram plot of the British coal mine disasters','fontsize',15)
xlabel('Year','fontsize',15)
ylabel('Frequency','fontsize',15)

for d = d_vek
    disp(['----- ANALYSIS FOR ',num2str(d-1),' BREAKPOINTS -----'])
    
    % We now try to sample from the distribution.
    
    %rho = 0.05;% + 0.10*(d-2);
    
    rho = rho_vek(d == d_vek);
    
    if d <= 4
        k = 1;bajs = d-1;
    else
        k = 2;bajs = 2;
    end
    
    lambda_vek = zeros(N,d);
    theta_vek = zeros(N,1);
    t_vek = zeros(N,d-1);

    t = sort([t_1 disa(randi([2,length(disa)-1],1,d-1)) t_end]);
    lambda = ones(1,d);
    
    func_t = @(t,lambda,n) prod(t)*exp(-sum(lambda.*...
        diff(t)))*prod(lambda.^n); 
    
    for i = 1:N

        n = zeros(1,d);
        for j = 1:d
            n(j) = length(disa(t(j) < disa & t(j+1) >= disa));
        end

        theta = gamrnd(2*(d + 1),1/(beta + sum(lambda)));
        lambda = gamrnd(2 + n,1./(theta + diff(t)));

        t_star = t;
        n_star = n;
        for j = 2:d
            
            
            R = rho*(t(j+1) - t(j-1));
            t_star(j) = t(j) - R + 2*R*rand;
            
            for l = 1:d
                n_star(l) = length(disa(t_star(l) < ...
                    disa & t_star(l+1) >= disa));
            end
            
            if t_star(j) > t(j-1) && t_star(j) < t(j+1)
           
                alpha = func_t(t_star,lambda,n_star)/func_t(t,lambda,n);
                helper = min(1,alpha);
                if rand < helper
                    t(j) = t_star(j);
                    n(j) = n_star(j);
                else
                    t_star(j) = t(j);
                    n_star(j) = n(j);
                end
            else
                t_star(j) = t(j);
                n_star(j) = n(j);
            end
        end

        lambda_vek(i,:) = lambda;
        theta_vek(i) = theta;
        t_vek(i,:) = t(2:end-1);
    end
   
    figure, hold on
    for j = 1:size(t_vek,2)
        %subplot(k,bajs,j),
        [temp1,temp2] = hist(t_vek(burn_in:end,j),20);
        bar(temp2,temp1/sum(temp1),color{j+1},'barwidth',0.3), hold on
        title(['Histogram plot of the estimated time for breakpoint '...
            num2str(j)],'fontsize',15)
        xlabel('Year','fontsize',15)
        ylabel('Percentage','fontsize',15)
    end
    
    figure
    for j = 1:size(lambda_vek,2)
        subplot(k,bajs+1,j), hist(lambda_vek(burn_in:end,j),20)
        title(['\lambda for i = ',num2str(j)]...
            ,'fontsize',15)
        xlabel('\lambda','fontsize',15)
        ylabel('Frequency','fontsize',15)
    end
    
    figure, hist(theta_vek(burn_in:end),20)
    title(['Plot of  \theta for d = ',num2str(d)],'fontsize',15)
    xlabel('\theta','fontsize',15)
    ylabel('Frequency','fontsize',15)
    
    figure
    for j = 1:size(t_vek,2)
        subplot(d-1,1,j), plot(t_vek(burn_in:end,j))
        title(['Plot of breakpoint ',num2str(j),' for d = ',num2str(d)],...
            'fontsize',25)
        xlabel('Sample','fontsize',25)
        ylabel('Value','fontsize',25)
    end

    for j = 1:size(t_vek,2)
        disp(['The year for breakpoint #',num2str(j),' is ',...
            num2str(mean(t_vek(burn_in:end,j)))])
    end
    
    disp(' ')
    
    disp(['Theta is ',num2str(mean(theta_vek(burn_in:end)))])
    
    disp(' ')
    
    for j = 1:size(lambda_vek,2)
        disp(['The intensity for interval #',num2str(j)...
            ,' is ',num2str(mean(lambda_vek(burn_in:end,j)))])
    end
    disp(' ')
    
    
    %if time, plot the distrubtions using the estimated parameter
    meant = [t_1 mean(t_vek(burn_in:end,:)) t_end];
    figure, bar(years_plot,disaster_plot), hold on
    for j = 2:d
        year = meant(j);
        yax = max(disaster_plot);
        plotters = 0:yax;
        plot(ones(1,yax+1)*year,plotters,color{j},'linewidth',1.5)
    end
    title(['The disasters divided into intervals for d = ',num2str(d)],...
        'fontsize',15)
    xlabel('Year','fontsize',15)
    ylabel('Frequency','fontsize',15)
    
    time = [];prutt = [];
    for j = 1:size(lambda_vek,2)
        time = [time meant(j:j+1)];
        prutt = [prutt mean(lambda_vek(burn_in:end,j))*ones(1,2)];
    end
    figure, plot(time,prutt,'--o')
    axis([min(meant) max(meant) 0 max(lambda_vek(:,1))])
    title(['Plot of  \lambda for d = ',num2str(d)],'fontsize',15)
    xlabel('Years','fontsize',15)
    ylabel('\lambda','fontsize',15)
end

