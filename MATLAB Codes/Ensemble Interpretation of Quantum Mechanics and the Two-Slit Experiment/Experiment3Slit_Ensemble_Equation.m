% Goal: Use the Ensemble Model to predict the outcome of the 3-slit
% experiment

% (0) Setting up interval
n = 2^10 + 1;                   % Number of points
X = linspace(-16,16, n);        % Partition
h = 2^(5 - 10);                 % Step size
assert(X(2) - X(1) == h)


% (1) Get the initial distribution of particle positions
s = 1;          % s = distance of means from origin
a = .1;         % a = standard deviation of means

W_0 = arrayfun(@(x)  InitialDistribution3Slit(x,s,a), X); 
max_x = 16;
max_y = max(W_0);
plot(X,W_0)
title("Initial Distribution 3-Slit Experiment")
xlim([-max_x max_x])
ylim([0 1])

% (2) Setting up the Model Parameters
alpha = 1 / (8*pi);
beta = 1 / (4*a)^2;

% (3) Numerically approximate solution
Z_end = 1 / pi;         % Distance interval of solution
k = 0.5*h^2;            % Step size in time
n_steps = round(Z_end / k);
W = W_0;
P1 = round(s/h);                % Number of steps to right or left to travel 1 unit

% (3.1) March forward in space
for step = [1:n_steps]
    for i = [1:n]
        % Update ith entry
        W_New(i) = W(i);  
        A=0;B=0;C=0;D=0;E=0;

        if i - P1 >= 1
            A = W(i - P1);
        end
        if i - 1 >= 1
            B = W(i - 1);
        end
        C = W(i);
        if i + 1 <= n
            D = W(i+1);
        end
        if i + P1 <= n
            E = W(i+P1);
        end
        derivative_Z = (alpha/h^2)*(B + D - 2*C) + beta*(A + E - 2*C);
        W_New(i) = W_New(i) + k * derivative_Z;
    end
    W_old = W;
    W = W_New;
    % (3.2) Plot our new pattern
    plot(X,W)
    title( sprintf("Ensemble Equation for the Double-Slit (z=%.3f)", round(step*k,3)) )
    xlim([-10 10])
    ylim([0 1])
    xlabel("x-position of particle")
    ylabel("Probability")
    Pause_Time = 30 / n_steps;         % (Total time) / (Number of time steps)
    pause(Pause_Time)
end




function [p] = InitialDistribution3Slit(x, s, a)
%InitialDistribution3Slit Return  initial distribution of particles
%   x = position along x-axis
%   s = distance between slits
%   a = initial standard deviation of Gaussian means 
const = 1 / (3*a*sqrt(2*pi));
E1 = exp(-1/2*  ((x+s) / a)^2 );
E2 = exp(-1/2*  ((x) / a)^2 );
E3 = exp(-1/2*  ((x-s) / a)^2 );
p = const * (E1 + E2 + E3);

if isnan(p)     % Fix case where p is really small!
    p = 0;
end

end
