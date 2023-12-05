% Problem 3 %
% ========= %

a = 3; % second last digit of roll number
b = 5; % last digit of roll number
L =1.4 ; % Length of both n and p region in micrometer
epsilon_Si = 11.8; % Dielectric constant of Silicon
epsilon_0 = 8.85 * 10^(-14); % Permittivity of free space in F/cm
n_i = 1.5*10^10; % Intrinsic carrier concentration per cm^3
E_g = 1.12; % Band gap energy in eV
T = 300; % Temperature in K
q = 1.6*10^(-19); % Charge of an electron
V_t = 0.026; % Thermal voltage in volts
V_bi = 0.06 + (a*b)/300;
N_d1 = 10^14; % Doping concentration of n region
N_d2 = N_d1*exp(V_bi/V_t);% Doping concentration of n-plus region

%%% T = 300 K %%%

h = 0.01; % Step size in mirometer
X = -L:h:L; % X-axis

% Potential matrix
V = zeros(length(X),1);

% Charge Density matrix assuming complete ionization
rho_v = zeros(length(X),1);

% Finite Difference matrix
A = zeros(length(X),length(X));
A(1,1) = 1;
A(1,2) = -1;
A(length(X),length(X)) = 1;
A(length(X),length(X)-1) = -1;
for i = 2:length(X)-1
    A(i,i-1) = 1;
    A(i,i) = -2;
    A(i,i+1) = 1;
end

% Matrix for dB/dV
B_1 = zeros(length(X),length(X));

for k = 1:500

    V(1,1) = V_t*log(N_d2/n_i); % Boundary Condition for n-plus region end
    V(length(X),1) = V_t*log(N_d1/n_i); % Boundary Condition for n region end
    
    p = n_i * exp(-V/V_t);
    n = n_i * exp(V/V_t);

    % Charge Density Matrix
    rho_v(1,1) = 0; % Boundary Condition for n-plus region end
    rho_v(2:round(L/h)+1,1) = q*(p(2:round(L/h)+1,1) - n(2:round(L/h)+1,1) + N_d2);
    rho_v(round(L/h)+2:length(X)-1,1) = q*(p(round(L/h)+2:length(X)-1,1) - n(round(L/h)+2:length(X)-1,1) + N_d1);
    rho_v(length(X),1) = 0; % Boundary Condition for n region end
    
    % Newton Raphson using Jacobian
    B = -(h^2*10^(-8)*(epsilon_Si*epsilon_0)^(-1))*rho_v;
    f = A * V - B;
    
    for i = 1:length(X)
        B_1(i,i) = -q*(h^2*10^(-8)*(epsilon_Si*epsilon_0)^(-1))*(V_t)^(-1)*n_i*(-exp(-V(i,1)/V_t)-exp(V(i,1)/V_t));
    end

    J = A - B_1;

    % L-U Decomposition
    LW = eye(length(X));
    UP = J;
    for i = 2:length(X)
        for j = i:length(X)
            LW(j,i-1) = UP(j,i-1)/UP(i-1,i-1);
            UP(j,:) = UP(j,:) - (UP(j,i-1)/UP(i-1,i-1))*UP(i-1,:);
        end
    end

    y = zeros(length(X),1);
    for i = 1:length(X)
        temp = 0;
        for j = 1:length(X)
            if(i~=j)
                temp = temp + LW(i,j)*y(j,1);
            end
        end
        y(i,1) = f(i,1) - temp;
    end
    
    x = zeros(length(X),1);
    for i = length(X):-1:1
        temp = 0;
        for j = 1:length(X)
            if(i~=j)
                temp = temp + UP(i,j)*x(j,1);
            end
        end
        x(i,1) = (y(i,1) - temp)/UP(i,i);
    end

    V = V - x;

end

% Charge Density in C/cm^3
p = n_i * exp(-V/V_t);
n = n_i * exp(V/V_t);
rho_v(1:round(L/h)+1,1) = q*(p(1:round(L/h)+1,1) - n(1:round(L/h)+1,1) + N_d2);
rho_v(round(L/h)+2:length(X),1) = q*(p(round(L/h)+2:length(X),1) - n(round(L/h)+2:length(X),1) + N_d1);
rho_v = rho_v';

% Potential in volts
V = V';

% Electric Field Intensity in V/cm
E = zeros(1,length(X));
for i = 2:length(X)
    E(1,i) = E(1,i-1) + h * 0.5 * (rho_v(1,i-1) + rho_v(1,i));
end
E = 10^(-4) * (epsilon_Si * epsilon_0)^(-1) * E;

% Potential, Electric Field, Charge Concentration Profile

figure;
plot(X* 10^(-4),V);
grid;
xlabel('x ( cm )');
ylabel('V ( volt )');
title('Voltage Profile');

fprintf("\nAt T = 300 K, Built-in potential = %f",V(1,1)-V(1,length(X)));

figure;
plot(X* 10^(-4),E);
grid;
xlabel('x ( cm )');
ylabel('E ( V/cm )');
title('Electric Field Profile');

figure;
plot(X* 10^(-4),rho_v/q);
grid;
xlabel('x ( cm )');
ylabel('Charge Concentration ( /cm^3 )');
title('Charge Concentration Profile');

% Electron and Hole concentration

p = n_i * exp(-V/V_t);
n = n_i * exp(V/V_t);

figure;
semilogy(X* 10^(-4),p,'Displayname','hole');
hold on;
semilogy(X* 10^(-4),n,'Displayname','electron');
xlabel('x ( cm )');
ylabel('Carrier Concentration ( /cm^3 )');
title('Electron and Hole Concentration');
legend;
grid;
hold off;

% Energy Band Diagram

E_F = zeros(1,length(X));
E_i = E_F - V;
E_C = E_i + E_g/2;
E_V = E_i - E_g/2;

figure;
plot(X,E_C,'Displayname','E_C');
hold on
plot(X,E_i,"--",'Displayname','E_{mid}');
plot(X,E_F,"k",'Displayname','E_F');
plot(X,E_V,"r",'Displayname','E_V');
ebd = gca;
ebd.XAxis.Visible = 'off';
ylabel("(E - E_F) in eV");
title("Energy Band Diagram");
legend;
grid;
hold off;

%%% T = 300+(1.5*a) K %%%

T1 = 300+(1.5*a);
V_t1 = 0.026*(T1/T);

% Potential matrix
V1 = zeros(length(X),1);

% Charge Density matrix assuming complete ionization
rho_v = zeros(length(X),1);

% Matrix for dB/dV
B_1 = zeros(length(X),length(X));

for k = 1:1000

    V1(1,1) = V_t1*log(N_d2/n_i); % Boundary Condition for n-plus region end
    V1(length(X),1) = V_t1*log(N_d1/n_i); % Boundary Condition for n region end
    
    p = n_i * exp(-V1/V_t1);
    n = n_i * exp(V1/V_t1);

    % Charge Density Matrix
    rho_v(1,1) = 0; % Boundary Condition for n-plus region end
    rho_v(2:round(L/h)+1,1) = q*(p(2:round(L/h)+1,1) - n(2:round(L/h)+1,1) + N_d2);
    rho_v(round(L/h)+2:length(X)-1,1) = q*(p(round(L/h)+2:length(X)-1,1) - n(round(L/h)+2:length(X)-1,1) + N_d1);
    rho_v(length(X),1) = 0; % Boundary Condition for n region end
    
    % Newton Raphson using Jacobian
    B = -(h^2*10^(-8)*(epsilon_Si*epsilon_0)^(-1))*rho_v;
    f = A * V1 - B;
    
    for i = 1:length(X)
        B_1(i,i) = -q*(h^2*10^(-8)*(epsilon_Si*epsilon_0)^(-1))*(V_t1)^(-1)*n_i*(-exp(-V1(i,1)/V_t1)-exp(V1(i,1)/V_t1));
    end

    J = A - B_1;

    % L-U Decomposition
    LW = eye(length(X));
    UP = J;
    for i = 2:length(X)
        for j = i:length(X)
            LW(j,i-1) = UP(j,i-1)/UP(i-1,i-1);
            UP(j,:) = UP(j,:) - (UP(j,i-1)/UP(i-1,i-1))*UP(i-1,:);
        end
    end

    y = zeros(length(X),1);
    for i = 1:length(X)
        temp = 0;
        for j = 1:length(X)
            if(i~=j)
                temp = temp + LW(i,j)*y(j,1);
            end
        end
        y(i,1) = f(i,1) - temp;
    end
    
    x = zeros(length(X),1);
    for i = length(X):-1:1
        temp = 0;
        for j = 1:length(X)
            if(i~=j)
                temp = temp + UP(i,j)*x(j,1);
            end
        end
        x(i,1) = (y(i,1) - temp)/UP(i,i);
    end

    V1 = V1- x;

end

V1 = V1';


figure;
plot(X* 10^(-4),V,'Displayname','T = 300 K');
hold on;
plot(X* 10^(-4),V1,'Displayname','T = 304.5 K');
xlabel('x ( cm )');
ylabel('V ( volt )');
title('Potential');
legend;
grid;
hold off;

fprintf("\nAt T = 304.5, Built-in potential = %f",V1(1,1)-V1(1,length(X)));
