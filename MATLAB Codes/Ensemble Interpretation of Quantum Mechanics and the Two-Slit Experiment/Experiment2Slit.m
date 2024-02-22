% Goal: Use the Schrodinger Model to predict outcome of the double-slit
% experiment

Schrodinger2Slit(1, 0, .1)
pause(10);
% (1) Set up initial distribution 
s = 1;          % s = distance of Gaussian means from origin
a = .1;         % a = standard deviation of means

n = 10000;
L = 10;                          % Length of Interval 
X = linspace(-L,L, n);
Y = arrayfun(@(x)  Schrodinger2Slit(x,0,a), X); 
max_y = max(Y);
plot(X,Y)
title("Initial Distribution 2-Slit Experiment")

% (2) Create animation show how interference pattern emerges
n = 1000;                         % Number of points
X = linspace(-L,L, n);      
n_steps = 300;                  % Number of time steps
end_time = 1/pi;

for t = linspace(0, end_time , n_steps)  
    % 1) Compute function
    Y = arrayfun(@(x)  Schrodinger2Slit(x,t,a), X); 

    % 2) Plot function
    plot(X,Y)                                           % 2) Plot
    title( sprintf("Distribution of Positions (t=%.3f)", round(t,3)) )
    xlim([-L, L])
    ylim([0, .5])
    xlabel("x-position of particle")
    ylabel("Probability")
    Pause_Time = 30 / n_steps;         % (Total time) / (Number of time steps)
    pause(Pause_Time)

    % 3) Check function
    % int = integral( @(x)  Schrodinger2Slit(x,t,a), 0, L);
    % assert( abs(int - 1/2) < 1e-3);
end




function [p] = Schrodinger2Slit(x, t, a)
%Schrodinger1Slit Solves the Schrodinger Equation for 1-slit
%   x = position along x-axis
%   t = time
%   a = initial standard deviation of Gaussian means 
N_0 = a.^4 + t.^2;
N_1  = exp( - (a.^2 + a.^2 * x.^2) / N_0);
N_2 = cosh( (2 * a.^2 * x) / N_0 ) + cos(2*t*x/N_0);
N = a .* N_1 .* N_2;
D = sqrt(pi) * sqrt(N_0) * (exp(1).^(-1/a.^2) + 1);
p = N / D;

if isnan(p)     % Fix case where p is really small!
    p = 0;
end

end

