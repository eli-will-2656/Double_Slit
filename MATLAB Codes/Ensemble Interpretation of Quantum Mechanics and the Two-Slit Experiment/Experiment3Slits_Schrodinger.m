% Goal: Predict the interference pattern shown in the 3-slit experiment


% (1.0) Select solution mesh 
% Don't make n too small or else erase initial condition
n = 100;                            % Number of points monitoring PDE at   
t_end = 1/pi;
x = linspace(-5,5, n);             % Interval monitoring interference pattern
t = linspace(0, t_end, round(t_end*10000));


% (1.1) Solve the equation
m = 0;
sol = pdepe(m,@pdefun,@pdeic,@pdebc,x,t);

% (1.2) Grab solutions for each function
u1 = sol(:,:,1);
u2 = sol(:,:,2);


surf(x,t,u1)
title('u_1(x,t)')
xlabel('Distance x')
ylabel('Time t')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) Set up the PDE
function [c,f,s] = pdefun(x,t,u,dudx)
c = [1; 1];
f = [1; -1] .* dudx;
s = [0; 0];
end

% (2) Set up initial condition
function u0 = pdeic(x)
u0 = [sqrt(InitialDistribution3Slit(x, 1, .1)); 0];
end

% (3) Set up boundary conditions
function [pl,ql,pr,qr] = pdebc(xl,ul,xr,ur,t)
pl = [0; ul(2)];
ql = [1; 0];        % Don't care about flux
pr = [ur(1); 0];
qr = [0; 1];        % Don't care about flux
end



% (4) Initial distribution
function [p] = InitialDistribution3Slit(x, s, a)
%InitialDistribution3Slit Return  initial distribution of particles
%   x = position along x-axis
%   s = distance between slits
%   a = initial standard deviation of Gaussian means 
const = 1 / (3*a*sqrt(2*pi));
E1 = exp(-1/2*  ((x+s) / a)^2 );
E2 = exp(-1/2*  (x/a)^2 );
E3 = exp(-1/2*  ((x-s)/a)^2 );
if isnan(E1)    
    E1 = 0;
end
if isnan(E2)   
    E2 = 0;
end
if isnan(E3)     
    E3 = 0;
end
p = const * (E1 + E2 + E3);

if isnan(p)     % Fix case where p is really small!
    p = 0;
end
end