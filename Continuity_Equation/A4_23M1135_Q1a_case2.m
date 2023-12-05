% Problem 1_A_Case_2 %
% ================== %

D = 0.1; % Diffusion coefficient in cm^2/s
n_0 = 10^12; % Concentration of particle at position A in cm^(-3)
taw = 10^(-7); % Minority carrier life time in sec 
k = 10^3; % Proportionality constant for flux cm/s
h = 0.000001; % Step size in cm
L = 10*10^(-4); % Length of the sample in cm
X = 0:h:L; % X-axis

A = zeros(length(X),length(X));
A(1,1) = 1; % Boundary Condition at x = 0
for i = 2:length(X)-1
    A(i,i-1) = 1/h^2;
    A(i,i) = -(2/h^2+(1/(D*taw)));
    A(i,i+1) = 1/h^2;
end
% Boundary Condition at x = 10 miron
A(length(X),length(X)-1) = -(1/h);
A(length(X),length(X)) = (1/h)+(k/D);

B = zeros(length(X),1);
B(1,1) = n_0;% Boundary Condition at x = 0
B(length(X),1) = 0; % Boundary Condition at x = 10 miron

n = A\B;
n = n';

% Particle Profile using Analytical Method

L_d = sqrt(D*taw); 
C_1 = n_0;
n_1 = C_1*exp(-X/L_d);

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

fprintf("\nParticle flux at B = %d cm^(-2)/s",k*n(1,length(X)));

figure;
semilogy(X,n);
xlabel('x ( cm )');
ylabel('Concentration ( /cm^3 )');
title('Particle Profile ( log scale )');
grid on;