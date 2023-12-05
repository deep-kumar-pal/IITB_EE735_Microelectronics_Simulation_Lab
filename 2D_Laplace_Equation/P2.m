% 2 %
% = %

L1 = 500; % length of the plate P1
d1 = 5; % distance between the plates P1 and P2/P3
d2 = [30 60 90]; % variable distance between plates P2 and P3
V3 = [0 5 -5]; % variable potential of plate P3
V = 5; % magnitude of potential of plate P1/P2
epsilon_0 = 8.854 * 10^(-12); % permittivity of free space (F/m)

% Implementing the experimental setup
m = (2 * L1) + 1; % No. of grids along x-axis
n = (3 * d1) + 1; % No. of grids along y-axis

% Some arbitary matrices
W = zeros(3,1);
C = zeros(3,3);
A = zeros(3,3);

for k=1:3 % no. of cases of distance
for l=1:3 % no. of cases of potentials

U = zeros(m,n); % Potential or Voltage matrix

% Solving Poisson's equation using Iteration method
p = 3000; % No. of iterations

for z = 1:p
        
    for i=2:m-1

        for j=2:n-1      

                % Potential of plate P1
                U(round(m/2) - floor(L1/2):round(m/2) + ceil(L1/2) - 1,round(n/2) - floor(d1/2) - 1) = -V;

                % Potential of plate P2
                U(round(m/2) - floor(L1/2):round(m/2) - round(d2(k)/2)-1,round(n/2) + ceil(d1/2)) = V;

                % Potential of plate P3
                U(round(m/2) + round(d2(k)/2):round(m/2) + ceil(L1/2) - 1,round(n/2) + ceil(d1/2)) = V3(l);
                
                % Implementing finite difference condition
                U(i,j)=0.25*(U(i+1,j)+U(i-1,j)+U(i,j+1)+U(i,j-1));
        end

     end
        
end
U = U';

% Calculating Electric Field
[Ex,Ey]=gradient(U);
Ex = -Ex;
Ey = -Ey;

% Electric field Magnitude
E = sqrt(Ex.^2+Ey.^2);

% Electrostatic Potential and Equipotential Surfaces plot
figure
pcolor(U);
hold on;
shading interp;
colorbar;
xlabel('x-axis in nanometers','fontsize',14);
ylabel('y-axis in nanometers','fontsize',14);
title(['ELECTROSTATIC POTENTIAL for d2 = ',num2str(d2(k)),' nm and V3 = ',num2str(V3(l)),' V']);

figure
contour(U);
colorbar;
xlabel('x-axis in nanometers','fontsize',14);
ylabel('y-axis in nanometers','fontsize',14);
title(['EQUIPOTENTIAL SURFACES for d2 = ',num2str(d2(k)),' nm and V3 = ',num2str(V3(l)),' V']);

% 2D Electric Field Profile plot
figure
pcolor(E);
hold on;
shading interp;
colorbar;
xlabel('x-axis in nanometers','fontsize',14);
ylabel('y-axis in nanometers','fontsize',14);
title(['2D ELECTRIC FIELD PROFILE for d2 = ',num2str(d2(k)),' nm and V3 = ',num2str(V3(l)),' V']);

% Energy per unit width
W(l,1) = 0.5*epsilon_0*sum(E.^2,"all");

V1 = -V; % potential of plate P1
V2 = V; % potential of plate P2
V12 = V1-V2; % potential difference between plates P1 and P2
V23 = V2-V3(l); % potential difference between plates P2 and P3
V31 = V3(l)-V1; % potential difference between plates P3 and P1

A(l,1) = 0.5*V12^2;
A(l,2) = 0.5*V23^2;
A(l,3) = 0.5*V31^2;

end

C(:,k) = A\W;

fprintf("\nFor d2 = %d nm :",d2(k));
fprintf("\nThe value of capacitance (per unit width), C_12 = %d F/m",C(1,k));
fprintf("\nThe value of capacitance (per unit width), C_23 = %d F/m",C(2,k));
fprintf("\nThe value of capacitance (per unit width), C_31 = %d F/m",C(3,k));

end