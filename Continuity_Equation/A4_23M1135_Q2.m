% PROBLEM 2 %
% ========= %

L = 100*10^(-4); % Length of the sample in cm
D = 10^(-4); % Diffusion coefficient in cm^2/s
n_0 = 2000; % Density at x = 0
t_final = 15 * 10^(-3); % Final observation at 15 millisec
t_0 = 10^(-3); % First observation at 1 millisec
h = 0.0001; % Step size along X-axis in cm
X = 0:h:L; % X-Axis

% Implicit or Backward Euler Method (C = 0.6)

C = 0.6; % Constant for Euler method
del_t = (C*h^2)/D; % Time spacing in sec
T = 0:del_t:t_final; % Time-axis

n = zeros(length(X),length(T)); % Concentration matrix
n(1,:) = n_0; % Forcing Constant particle at x = 0

for k = 1:4000
    for t = 2:length(T)
        for i = 2:length(X)-1
            % Backward Euler Condition
            n(i,t) = (C*(n(i-1,t)+n(i+1,t))+n(i,t-1))/(2*C+1);
        end
    end
end

figure;
plot(X,n(:,floor(t_0/del_t+1)),'Displayname','t = t_0');
hold on;
plot(X,n(:,floor(2*t_0/del_t+1)),'Displayname','t = 2t_0');
plot(X,n(:,floor(3*t_0/del_t+1)),'Displayname','t = 3t_0');
plot(X,n(:,floor(4*t_0/del_t+1)),'Displayname','t = 4t_0');
plot(X,n(:,floor(5*t_0/del_t+1)),'Displayname','t = 5t_0');
plot(X,n(:,floor(6*t_0/del_t+1)),'Displayname','t = 6t_0');
plot(X,n(:,floor(7*t_0/del_t+1)),'Displayname','t = 7t_0');
plot(X,n(:,floor(8*t_0/del_t+1)),'Displayname','t = 8t_0');
plot(X,n(:,floor(9*t_0/del_t+1)),'Displayname','t = 9t_0');
plot(X,n(:,floor(10*t_0/del_t+1)),'Displayname','t = 10t_0');
plot(X,n(:,floor(11*t_0/del_t+1)),'Displayname','t = 11t_0');
plot(X,n(:,floor(12*t_0/del_t+1)),'Displayname','t = 12t_0');
plot(X,n(:,floor(13*t_0/del_t+1)),'Displayname','t = 13t_0');
plot(X,n(:,floor(14*t_0/del_t+1)),'Displayname','t = 14t_0');
plot(X,n(:,floor(15*t_0/del_t+1)),'Displayname','t = 15t_0');
title("Evolution of Particle Density (Numerical)");
xlabel("x ( cm )");
ylabel("Concentration ( /cm^3 )");
legend;
grid;
hold off;

% Analytical Method

n_1 = zeros(length(X),length(T)-1); % Concentration matrix
for k = 1:length(T)-1
    % Profile Function w.r.t space and time
    n_1(:,k) = n_0*erfc(X/(2*sqrt(D*T(1,k+1))));
end

figure;
plot(X,n_1(:,floor(t_0/del_t)),'Displayname','t = t_0');
hold on;
plot(X,n_1(:,floor(2*t_0/del_t)),'Displayname','t = 2t_0');
plot(X,n_1(:,floor(3*t_0/del_t)),'Displayname','t = 3t_0');
plot(X,n_1(:,floor(4*t_0/del_t)),'Displayname','t = 4t_0');
plot(X,n_1(:,floor(5*t_0/del_t)),'Displayname','t = 5t_0');
plot(X,n_1(:,floor(6*t_0/del_t)),'Displayname','t = 6t_0');
plot(X,n_1(:,floor(7*t_0/del_t)),'Displayname','t = 7t_0');
plot(X,n_1(:,floor(8*t_0/del_t)),'Displayname','t = 8t_0');
plot(X,n_1(:,floor(9*t_0/del_t)),'Displayname','t = 9t_0');
plot(X,n_1(:,floor(10*t_0/del_t)),'Displayname','t = 10t_0');
plot(X,n_1(:,floor(11*t_0/del_t)),'Displayname','t = 11t_0');
plot(X,n_1(:,floor(12*t_0/del_t)),'Displayname','t = 12t_0');
plot(X,n_1(:,floor(13*t_0/del_t)),'Displayname','t = 13t_0');
plot(X,n_1(:,floor(14*t_0/del_t)),'Displayname','t = 14t_0');
plot(X,n_1(:,floor(15*t_0/del_t)),'Displayname','t = 15t_0');
title("Evolution of Particle Density (Analytical)");
xlabel("x ( cm )");
ylabel("Concentration ( /cm^3 )");
legend;
grid;
hold off;

% Analytical vs. Numerical

figure;
plot(X,n_1(:,floor(t_0/del_t)),'Displayname','Analytical ( t = t_0 )');
hold on;
plot(X,n(:,floor(t_0/del_t+1)),'Displayname','Numerical ( t = t_0 )');
plot(X,n_1(:,floor(2*t_0/del_t)),'Displayname','Analytical ( t = 2t_0 )');
plot(X,n(:,floor(2*t_0/del_t+1)),'Displayname','Numerical ( t = 2t_0 )');
plot(X,n_1(:,floor(3*t_0/del_t)),'Displayname','Analytical ( t = 3t_0 )');
plot(X,n(:,floor(3*t_0/del_t+1)),'Displayname','Numerical ( t = 3t_0 )');
title("Analytical vs. Numerical");
xlabel("x ( cm )");
ylabel("Concentration ( /cm^3 )");
legend;
grid;
hold off;