% ----- ERROR ANALYSIS -----

% ----- ANALYSIS OF HYPERPARAMETER -----

disp(' ----- ANALYSIS OF HYPERPARAMETER -----')
disp(' ')

d_vek = 3; 

beta_vek = [0 1 5 20];

% This parameter represents the number of intervals we have, i.e. if we 
% have d = 3 we have two breakpoints.

% From looking at the plots for different beta and rho, we can easily see
% that the MCMC algorithm is extremely robust for beta, wheras we can't
% vary rho that much.

% Tuning parameter

burn_in = 1e4;

rho = 0.14;

% We define the burn in to be 5000

t_1 = 1851;t_end = 1963;N = 2.5e4 + burn_in;

% We then load the coal mine disasters

load coal_mine.mat

disa = coal_mine';

% We then plot the disasters 

runs = 100;

% We define the function

func_t = @(t,lambda,n) prod(t)*exp(-sum(lambda.*...
    diff(t)))*prod(lambda.^n); 

for d = d_vek
    schtek = [];
    for beta = beta_vek
        saver = 1;
        mat = zeros(runs,2*d);
        for run = 1:runs
            
            % We now try to sample from the distribution.

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
            mat(saver,:) = [mean(theta_vek(burn_in:end)) ...
                mean(t_vek(burn_in:end,:)) mean(lambda_vek(burn_in:end,:))];
            saver = saver + 1;
        end
        normz = var(mat);
        schtek = [schtek;normz];
        disp(['beta = ',num2str(beta)])
        disp(['The variance of theta is ',num2str(normz(1))])
        disp(['The variance of t is ',num2str(normz(2:d))])
        disp(['The variance of lambda is ',num2str(normz(d+1:end))])
        disp(' ')
    end
    
end