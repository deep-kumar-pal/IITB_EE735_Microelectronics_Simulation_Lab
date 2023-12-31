% 1(A,B,C,D) %
% ========== %


% Dimensions of the Capacitor
L = 800; % Length of the plates
d = 8; % Distance between the plates
V = 1; % Potential of plate
epsilon_0 = 8.854 * 10^(-12); % Permittivity of free space (F/m)
epsilon_r = 1; % Relative permittivity of medium

% Implementing the experimental setup
m = (2 * L) + 1; % No. of grids along x-axis
n = (3 * d) + 1; % No. of grids along y-axis

% Initialising the potential matrix
U = zeros(m,n);

% Solving Poisson's equation using Iteration method

p = 3000; % No. of iterations

for z = 1:p
        
    for i=2:m-1

        for j=2:n-1      
                
                % Potential of lower plate
                U(round(m/2) - floor(L/2) - 1:round(m/2) + ceil(L/2),round(n/2) - floor(d/2) - 1) = -V;
                % Potential of upper plate
                U(round(m/2) - floor(L/2) - 1:round(m/2) + ceil(L/2),round(n/2) + ceil(d/2)) = V;
                % Implementation of FDM
                U(i,j) = 0.25*(U(i+1,j)+U(i-1,j)+U(i,j+1)+U(i,j-1));
        end

    end
        
end
U = U';


% Electric Field
[Ex,Ey]=gradient(U);
Ex = -Ex;
Ey = -Ey;

% Electric field Magnitude
E = sqrt(Ex.^2+Ey.^2);

% Energy per unit width (J/m)
W = 0.5*epsilon_r*epsilon_0*sum(E.^2,"all");

% Capacitance per unit width (F/m)
C = W/(2*V^2);

fprintf("\nThe value of capacitance (per unit width) is, C = %d F/m",C);

% Electrostatic Potential and Equipotential Surfaces
figure
contour(U);
colorbar;
xlabel('x-axis in nanometers','fontsize',14);
ylabel('y-axis in nanometers','fontsize',14);
title("EQUIPOTENTIAL SURFACES");

figure
pcolor(U);
hold on;
shading interp;
colorbar;
xlabel('x-axis in nanometers','fontsize',14);
ylabel('y-axis in nanometers','fontsize',14);
title("ELECTROSTATIC POTENTIAL");

% 2D Electric Field Profile
figure
pcolor(E);
hold on;
shading interp;
colorbar;
xlabel('x-axis in nanometers','fontsize',14);
ylabel('y-axis in nanometers','fontsize',14);
title("2D ELECTRIC FIELD PROFILE");

% Theoretical value of capacitance per unit width (F/m)
C_th = epsilon_r*epsilon_0*(L/d);
fprintf("\nThe value of theoretical capacitance (per unit width) is, C_th = %d F/m",C_th);