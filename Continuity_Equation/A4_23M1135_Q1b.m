% Problem 1_B %
% =========== %

D = 0.1; % Diffusion coefficient in /cm^2-s
L = 10*10^(-4); % Length of the sample in cm
p = (3+0.5*5)*10^(-4); % Position of Injection in cm
J_0 = 10^13; % Flux near point of injection in /cm^2-s
taw = 10^(-7); % Minority carrier lifetime in sec
h = 0.0000001; % Step size in cm
X = 0:h:L; % X-axis

% Particle Profile using Numerical Analysis

A = zeros(length(X),length(X)); % Finite Difference Matrix
% Boundary condition at A and B
A(1,1) = 1;
A(length(X),length(X)) = 1;
for i = 2:length(X)-1
    A(i,i-1) = 1/h^2;
    A(i,i) = -(2/h^2+(1/(D*taw)));
    A(i,i+1) = 1/h^2;
end
% Boundary condition at point of injection 
A(round(p/h+1),round(p/h+1)) = 2/h;
A(round(p/h+1),round(p/h)) = -1/h;
A(round(p/h+1),round(p/h+2)) = -1/h;

B = zeros(length(X),1);

% Boundary condition at A, B and point of injection
B(1,1) = 0;
B(length(X),1) = 0;
B(round(p/h+1),1) = 2*J_0/D;

% Particle Concentration matrix
n = (A\B);
n = n';

figure;
plot(X,n);
xlabel('x ( cm )');
ylabel('Concentration ( /cm^3 )');
title('Particle Profile : Numerical Analysis');
grid;

% Flux Profile
J = -D*gradient(n,X);
fprintf("\nParticle flux at A = %d cm^(-2)/s",J(1,1));
fprintf("\nParticle flux at B = %d cm^(-2)/s",J(1,length(X)));


% Particle Profile using Analytical Method

L_d = sqrt(D*taw);
n_1 = zeros(1,length(X));

% Particle profile from A to point of injection
C_1 = -(10^10/(exp(-5.5)+exp(5.5)));
C_2 = 10^10/(exp(-5.5)+exp(5.5));
n_1(1,1:round(p/h+1)) = C_1*exp(-X(1,1:round(p/h+1))/L_d)+C_2*exp(X(1,1:round(p/h+1))/L_d);

% Particle profile from point of injection to B
C_3 = 10^10/(1+exp(-9));
C_4 = -(10^10/(1+exp(9)));
n_1(1,round(p/h+2):length(X)) = C_3*exp(-(X(1,round(p/h+2):length(X))-p)/L_d)+C_4*exp((X(1,round(p/h+2):length(X))-p)/L_d);

figure;
plot(X,n,'Displayname','Numerical');
hold on;
plot(X,n_1,'Displayname','Analytical','LineStyle','--');
xlabel('x ( cm )');
ylabel('Concentration ( /cm^3 )');
title('Particle Profile : Numerical vs. Analytical');
legend;
grid;
hold off;