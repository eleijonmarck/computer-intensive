% ----- TASK 2--------
disp('----- TASK 2 -----')
disp(' ')

% We begin by importing the waves and get that

load atlantic.mat

% Nameing it as waves

waves = atlantic;

% We then plot the waves

plot((0:length(waves)-1)/length(waves)*100,waves,'--')
title('Plot of the waves','fontsize',25)
xlabel('Year','fontsize',25)
ylabel('Wave height','fontsize',25)

% We then get the estimates for the parameters by using est_gumbel

[beta,mu] = est_gumbel(waves);

% We then sample from the Gumbel distribution by defining the function

gumbrand = @(x,beta,mu) -beta*log(log(1./x)) + mu;

% Now we begin the bootstrapping

n = length(waves);
B = 5e2;
boots = zeros(B,2);
for i = 1:B
    X = gumbrand(rand(1,n),beta,mu);
    [bb,bm] = est_gumbel(X);
    boots(i,:) = [bb bm];
end
deltas = sort(bsxfun(@minus,boots,[beta mu]));
L = bsxfun(@minus,[beta mu],deltas(ceil((1 - 0.05/2)*B),:));
U = bsxfun(@minus,[beta mu],deltas(ceil(0.05/2*B),:));

disp('Bounds for beta and mu of the estimated Gumpel')
disp(' ')
disp('      L         U')
disp([L' U'])

% We now defined T to be

T = 3*14*100;

% We now calculate the maximum value which is 

max_wave = gumbrand(1 - 1/T,U(1),U(2));
disp(['Maximum wave height for a 100-year period is ',num2str(max_wave),...
    ' m'])
disp(['The confidence bound is therefore (',num2str(gumbrand(1 - 1/T,...
    beta,mu)),',',num2str(max_wave),')'])

gumpdf = @(x,beta,mu) 1/beta*exp(-((x - mu)/beta + exp(-(x - mu)/beta)));

figure, bar(hist(waves,14)/sum(hist(waves,14))), hold on
plot(0:0.01:max(waves),gumpdf(0:0.01:max(waves),U(1),U(2)),'r--'),
plot(0:0.01:max(waves),gumpdf(0:0.01:max(waves),beta,mu),'k--'),
plot(0:0.01:max(waves),gumpdf(0:0.01:max(waves),L(1),L(2)),'g--'),
title('Histogram plot of the sampled gumpel distribution', 'Fontsize', 25)
xlabel('Value', 'Fontsize', 25)
ylabel('Frequency', 'Fontsize', 25)
legend('Raw data','Upper Gumbel PDF', 'Lower Gumpel PDF')
