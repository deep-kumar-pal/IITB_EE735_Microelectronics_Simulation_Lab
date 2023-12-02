% 1(E) %
% ==== %

% Dimensions of the Capacitor
L = 800; % Length of the plates
d = 8; % Distance between the plates
V = 1; % Potential of plate
epsilon_0 = 8.854 * 10^(-12); % Permittivity of free space (F/m)
epsilon_r1 = 1; % Relative permittivity of first medium
epsilon_r2 = 5; % Relative permittivity of second medium

% Implementing the experimental setup
m = (2 * L) + 1; % No. of grids along x-axis
n = (3 * d) + 1; % No. of grids along y-axis

% Initialising the potential matrix
U = zeros(m,n); % Potential or Voltage matrix

% Solving Poisson's equation using Iteration method
p = 3000; % No. of iterations

for z = 1:p
        
    for i=2:m-1

        for j=2:round(n/2)      

                % Potential of lower plate
                U(round(m/2) - floor(L/2) - 1:round(m/2) + ceil(L/2)-1,round(n/2) - floor(d/2)) = -V;
                
                % Implementing finite difference condition
                U(i,j)=0.25*(U(i+1,j)+U(i-1,j)+U(i,j+1)+U(i,j-1));
        end

        for j=round(n/2)+1:n-1

                % Potential of upper plate
                U(round(m/2) - floor(L/2) - 1:round(m/2) + ceil(L/2)-1,round(n/2) + ceil(d/2) +1) = V;
                
                % Implementing finite difference condition
                U(i,j)=0.25*(U(i+1,j)+U(i-1,j)+U(i,j+1)+U(i,j-1));
        end

                % Implementing Boundary Condition
                U(i,round(n/2))=(epsilon_r1*U(i,round(n/2)+2)-epsilon_r1*U(i,round(n/2)+1)+epsilon_r2*U(i,round(n/2)-1))/epsilon_r2;
     
    end
        
end
U = U';

% Dielectric matrix
epsilon_mat=zeros(n,m);
epsilon_mat(1:round(n/2),:)=epsilon_r2;
epsilon_mat(round(n/2)+1:n,:)=epsilon_r1;
epsilon_mat(round(n/2) + ceil(d/2) +1,round(m/2) - floor(L/2) - 1:round(m/2) + ceil(L/2)-1)=0;
epsilon_mat(round(n/2) - floor(d/2),round(m/2) - floor(L/2) - 1:round(m/2) + ceil(L/2)-1)=0;

% Calculating Electric Field for medium 1
[Ex1,Ey1]=gradient(U(round(n/2)+1:n,:));
Ex1 = -Ex1;
Ey1 = -Ey1;

% Calculating Electric Field for medium 2
[Ex2,Ey2]=gradient(U(1:round(n/2),:));
Ex2 = -Ex2;
Ey2 = -Ey2;

% Combined Electric Field
Ex = zeros(n,m);
Ey = zeros(n,m);
Ex(1:round(n/2),:) = Ex2;
Ex(round(n/2)+1:n,:) = Ex1;
Ey(1:round(n/2),:) = Ey2;
Ey(round(n/2)+1:n,:) = Ey1;

% Electric field Magnitude
E = sqrt(Ex.^2+Ey.^2);

% Energy per unit width (J/m)
W = 0.5*epsilon_0*sum(epsilon_mat.*(E.^2),"all");

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
C_1 = epsilon_r1*epsilon_0*((2*L)/d);
C_2 = epsilon_r2*epsilon_0*((2*L)/d);
C_th = (C_1*C_2)/(C_1+C_2);
fprintf("\nThe value of theoretical capacitance (per unit width) is, C_th = %d F/m",C_th);