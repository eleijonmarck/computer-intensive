clc, clear all, close all

load coal_mine.mat
dist = coal_mine;
t_1 = 1812; t_end = 1962;
N=200;

d = 2; %number of intervals
bp = d-1; %number of breakpoints
vi = 1; %hyperprior guess

t = sort([t_1 dist(randi([2,length(dist)-1],d-1)) t_end]);



for k = 1:bp
    
    for i=1:N %Gibbs sampling
        
        %theta = gamrnd(2*(d+1),1/(vi+sum(lambda)));
        %lambda = gamrnd(sum(dist)+2,1./(theta+diff(t)));
    end
end
