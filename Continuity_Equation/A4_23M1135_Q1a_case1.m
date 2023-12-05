% Problem 1_A_Case_1 %
% ================== %

D = 0.1; % Diffusion coefficient in cm^2/s
n_0 = 10^12; % Concentration of particle at position A in cm^(-3)
taw = 10^(-7); % Minority carrier life time in sec
h = 0.000001; % Step size in cm
L = 10*10^(-4); % Length of the sample in cm
X = 0:h:L; % X-axis

% Particle Profile using Numerical Method

A = zeros(length(X),length(X)); % Finite Difference Matrix
A(1,1) = 1; % Boundary Condition at x = 0
A(length(X),length(X)) = 1; % Boundary Condition at x = 10 miron

for i = 2:length(X)-1
    A(i,i-1) = 1/h^2;
    A(i,i) = -(2/h^2+(1/(D*taw)));
    A(i,i+1) = 1/h^2;
end

B = zeros(length(X),1);
B(1,1) = n_0; % Boundary Condition at x = 0
B(length(X),1) = 0; % Boundary Condition at x = 10 miron

n = A\B;
n = n';

% Particle Profile using Analytical Method

L_d = sqrt(D*taw); 
C_1 = n_0/(1-exp(-20));
C_2 = n_0/(1-exp(20));
n_1 = C_1*exp(-X/L_d)+C_2*exp(X/L_d);

figure;
plot(X,n_1,'Displayname','Analytical');
hold on;
plot(X,n,'Displayname','Numerical','LineStyle','--');
xlabel('x ( cm )');
ylabel('Concentration ( /cm^3 )');
title('Particle Profile : Numerical vs. Analytical');
legend;
grid;
hold off;

% Flux Profile

J = -D*gradient(n,X); % Numerical
J_1 = (D/L_d)*(C_1*exp(-X/L_d)-C_2*exp(X/L_d)); % Analytical

figure;
plot(X,J_1,'Displayname','Analytical');
hold on;
plot(X,J,'Displayname','Numerical','LineStyle','--');
xlabel('x ( cm )');
ylabel('Flux ( /cm^2-s )');
title('Flux Profile : Numerical vs. Analytical');
legend;
grid;
hold off;