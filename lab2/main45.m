% ----- TASK 1 -----

clc, clear all, close all

disp('----- TASK 1 -----')
disp(' ')

% We're just to test whether we can sample from one of the densities. Thus
% we begin by defining two parameters. 

d_vek = 2:5; 

% This parameter represents the number of intervals we have, i.e. if we 
% have d = 3 we have two breakpoints.

beta = 1; % GOOD BETA IS 1

% Hyper parameter

rho = 1; % GOOD HYPER PARAMETER IS 0.05

% From looking at the plots for different beta and rho, we can easily see
% that the MCMC algorithm is extremely robust for beta, wheras we can't
% vary rho that much.

% Tuning parameter

burn_in = 5e3;

% We define the burn in to be 5000

t_1 = 1851;t_end = 1963;N = 5e4 + burn_in;

% We then load the coal mine disasters

load coal_mine.mat

disa = coal_mine';

% We then plot the disasters 

hist(disa,20)
title('Histogram plot of the British coal mine disasters','fontsize',15)
xlabel('Year','fontsize',15)
ylabel('Frequency','fontsize',15)

func_t = @(t,lambda,n) prod(diff(t))*exp(-sum(lambda)*diff(t))*prod(lambda.^n); %Anders magi

for d = d_vek
    disp(['----- ANALYSIS FOR ',num2str(d-1),' BREAKPOINTS -----'])
    % We now try to sample from the distribution.
    
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
    
    for i = 1:N

        n = zeros(1,d);
        for j = 1:d
            n(j) = length(disa(t(j) < disa & t(j+1) >= disa));
        end

        theta = gamrnd(2*(d + 1),1/(beta + sum(lambda)));
        lambda = gamrnd(2 + n,1./(theta + diff(t)));

        t_help = t;
        t_star = t;
        for j = 2:d
            R = rho*(t(j+1) - t(j-1));
            t_star(j) = t(j) - R + 2*R*rand;
            
            for l = 1:d
                n_star(l) = length(disa(t_star(l) < disa & t_star(l+1) >= disa));
            end
            
            if t_star(j) > t(j-1) && t_star(j) < t(j+1)
           
                
                helper = min(1,func_t(t_star,lambda,n_star)/func_t(t,lambda,n));
                if rand <= helper
                    t = t_star;
                    n = n_star;
                else
                    t_star = t;
                end
              
            else
                t_star = t;
            end
        end

        lambda_vek(i,:) = lambda;
        theta_vek(i) = theta;
        t_vek(i,:) = t(2:end-1);
    end
   
    figure
    for j = 1:size(t_vek,2)
        subplot(k,bajs,j),hist(t_vek(burn_in:end,j),25)
        title(['Histogram plot of the estimated time for breakpoint '...
            num2str(j)],'fontsize',15)
        xlabel('Year','fontsize',15)
        ylabel('Frequency','fontsize',15)
    end

    for j = 1:size(t_vek,2)
        disp(['The year for breakpoint #',num2str(j),' is ',...
            num2str(mean(t_vek(:,j)))])
    end
    
    disp(' ')
    
    for j = 1:size(lambda_vek,2)
        disp(['The intensity for interval #',num2str(j)...
            ,' is ',num2str(mean(lambda_vek(:,j)))])
    end
    disp(' ')
    
end