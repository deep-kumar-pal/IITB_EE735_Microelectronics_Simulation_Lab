% 3 %
% = %


% Dimensions of the Capacitor
L = 20:50:1470; % Length of the plates
d = 8; % Distance between the plates
V = 1; % Potential of plate
epsilon_0 = 8.854 * 10^(-12); % permittivity of free space (F/m)
epsilon_r = 1; % relative permittivity of medium

q = length(L);
C = zeros(1,q);
for k = 1:q
% Implementing the experimental setup
m = (3 * L(k)) + 1; % No. of grids along x-axis
n = (3 * d) + 1; % No. of grids along y-axis

% Initialising the potential matrix
U = zeros(m,n); % Potential or Voltage matrix

% Solving Poisson's equation using Iteration method
p = 2000; % No. of iterations

for z = 1:p
        
        for i=2:m-1
        for j=2:n-1      
                
                % Potential of upper plate
                U(round(m/2) - floor(L(k)/2) - 1:round(m/2) + ceil(L(k)/2),round(n/2) - floor(d/2) - 1) = -V;
                % Potential of lower plate
                U(round(m/2) - floor(L(k)/2) - 1:round(m/2) + ceil(L(k)/2),round(n/2) + ceil(d/2)) = V;
                
                U(i,j)=0.25*(U(i+1,j)+U(i-1,j)+U(i,j+1)+U(i,j-1));
        end
        end
        
end
U = U';

% Calculating Electric Field
[Ex,Ey]=gradient(U);
Ex = -Ex;
Ey = -Ey;

% Calculation of Capacitance
E = sqrt(Ex.^2+Ey.^2); % Electric field Magnitude
W = 0.5*epsilon_r*epsilon_0*sum(E.^2,"all"); % Energy per unit width
C(1,k) = W/(2*V^2); % Capacitance per unit width (F/m)

fprintf("\nSimulation instance %d of %d completed.",k,q);

end

% C vs L
plot(L,C);
xlabel("Length of the plates (L in nm)");
ylabel("Capacitance per unit width (C in F/m)");
title("C vs L");

% C_p vs L
C_th = epsilon_r*epsilon_0*(L/d);
C_p = C - C_th;
plot(L,C_p);
xlabel("Length of the plates (L in nm)");
ylabel("Parasitic capacitance per unit width (C_p in F/m)");
title("C_p vs L");

% Relative Error
relative_error = C_p./C;
plot(L,relative_error);
xlabel("Length of the plates (L in nm)");
ylabel("C_p/C");
title("Relative Error vs L");
