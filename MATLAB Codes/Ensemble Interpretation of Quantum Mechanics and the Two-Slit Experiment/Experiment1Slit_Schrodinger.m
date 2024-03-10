% Goal: Recreate the Figures in Glenn's 2018 paper on the Ensemble Equation


% Section 1: Introduction 


% Section 2: Schrodinger Equation Model of the Double Slit Experiment

% (2.1) Setting up Initial Distribution of Particle Positions
s = 0;          % s = mean
a = 1;          % a = standard deviation

% Drawing initial distribution
n = 100;
X = linspace(-5,5, n);
Y = arrayfun(@(x)  Gaussian1Slit(x,s,a), X); 
plot(X,Y)
title("Initial Distribution 1-Slit Experiment")


% (2.2) Solve the Schrodinger Equation for this situation
n = 100;                        % Number of points
L = 10;                          % Length of Interval 
X = linspace(-L,L, n);      
n_steps = 100;                  % Number of time steps

for t = linspace(0, 10, n_steps)  
    % 1) Compute function
    Y = arrayfun(@(x)  Schrodinger1Slit(x,t,s,a), X);  

    % 2) Plot function
    plot(X,Y)                                           % 2) Plot
    title( sprintf("Distribution of Positions (t=%.3f)", round(t,3)) )
    xlim([-L, L])
    ylim([0, .5])
    xlabel("x-position of particle")
    ylabel("")
    Pause_Time = 5 / n_steps;         % (Total time) / (Number of time steps)
    pause(Pause_Time)

    % 3) Check function
    int = integral( @(x)  Schrodinger1Slit(x,t,s,a), 0, Inf);
    assert( abs(int - 1/2) < 1e-6);
end



function [p] = Gaussian1Slit(x, s, a)
%Gaussian1Slit Creates an initial distribution of particle positions 1-slit
%   x = position along x-axis
%   s = mean
%   a = standard deviation
D = (s-x).^2 / (2*a.^2);      % Quantifies distance from the mean
p = exp(-D) / (sqrt(2*pi) * a);

% Alternatively
% p = normpdf(x, s, a);
end

function [p] = Schrodinger1Slit(x, t, s, a)
%Schrodinger1Slit Solves the Schrodinger Equation for 1-slit
%   x = position along x-axis
%   t = time
%   s = mean
%   a = standard deviation
I = (sqrt(2)*a*(s-x)).^2 / (4*a.^4 + t.^2);        % Inside exponential
N = sqrt(2)*a * exp(-I);                           % Numerator
D = sqrt(pi)*sqrt(4*a^2 + t.^2);
p = N / D;
end
