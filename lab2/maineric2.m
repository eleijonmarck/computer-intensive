%HA2 main program

clc, clear all, close all

%Coal mine disasters-constructing a complex MCMC algorithm
%here we want to creat a HYBRID MCMC, i.e. Gibbs with MH


tot_dis = 191; %total number of disasters
bp = 1; %breakpoints
N = 200; %particles simulated 

burn_in = 2000;
X = 
X_k = [];

for i=1:bp

	burn_in = 2000;
	X = randn(1,burn_in+N); %initial xi distribution
	X_k = [X X_k]
	n_i = tau_i
	for k=1:N-1
		X_star =  proposal_kernel(1,N);%sample from proposal kernel	
		beta = func(X_star)*prop_kernel(
		alpha = min(1,beta);
%To do

%number of disasters at i:th interval
n_i = tau_i;

%generate porportional poisson processes for i:th interval
poissapdf = @(lambda,n_i) e^(-lambda)*lambda^n_i;

%density of disaster
disaster_density = poissapdf(n_i,lambda);

%all breakpoints are to be implemented with MH







%% Examples
% boring to do this way?
%Metropolis Hastings Matlab
rng('default');  % For reproducibility
delta = .5;
pdf = @(x) normpdf(x);
proppdf = @(x,y) unifpdf(y-x,-delta,delta);
proprnd = @(x) x + rand*2*delta - delta;
nsamples = 15000;
x = mhsample(1,nsamples,'pdf',pdf,'proprnd',proprnd,'symmetric',1);


%A gibb sampler example
burn_in = 1000;
M = N + burn_in;
X = zeros(1,M);
Y = X;
X(1) = 5;
Y(1) = 0.5;
for k = 1:(M - 1),
	x = binornd(n,Y(k));
	X(k + 1) = x;
	Y(k + 1) = betarnd(x + alpha,n - x + beta);
end

%A block technique

K = 50; % block size
n = N/K; % number of blocks
T = zeros(1,n);
for k = 1:n, % take means over n blocks
T(k) = mean(Y((burn_in + (k - 1)*K + 1):(burn_in + K*k)));
end
LB = tau - norminv(0.975)*std(T)/sqrt(n); % confidence bound
UB = tau + norminv(0.975)*std(T)/sqrt(n);

