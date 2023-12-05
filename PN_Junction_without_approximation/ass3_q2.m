% Problem 2 %
% ========= %

L_p1 = 0.5; % Length p-plus region in micrometer
L_n = 2; % Length n region in micrometer
L_n1 = 0.5; % Length n-plus region in micrometer
epsilon_Si = 11.8; % Dielectric constant of Silicon
epsilon_0 = 8.85 * 10^(-14); % Permittivity of free space in F/cm
n_i = 1.5 * 10^10; % Intrinsic carrier concentration per cm^3
E_g = 1.12; % Band gap energy in eV
T = 300; % Temperature in K
q = 1.6*10^(-19); % Charge of an electron
V_t = 0.026; % Thermal voltage in volts
N_a1 = 10^17; % Doping concentration of p-plus region
N_d = 10^15;% Doping concentration of n region
N_d1 = 10^17; % Doping concentration of n-plus region

h = 0.01; % Step size in mirometer
X = -L_p1:h:(L_n+L_n1); % X-axis

X2 = X;

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

for k = 1:1000

    V(1,1) = -V_t*log(N_a1/n_i); % Boundary Condition for n-plus region end
    V(length(X),1) = V_t*log(N_d1/n_i); % Boundary Condition for n region end
    
    p = n_i * exp(-V/V_t);
    n = n_i * exp(V/V_t);

    % Charge Density Matrix
    rho_v(1,1) = 0; % Boundary Condition for n-plus region end
    rho_v(2:round(L_p1/h)+1,1) = q*(p(2:round(L_p1/h)+1,1) - n(2:round(L_p1/h)+1,1) - N_a1);
    rho_v(round(L_p1/h)+2:round((L_p1+L_n)/h)+1,1) = q*(p(round(L_p1/h)+2:round((L_p1+L_n)/h)+1,1) - n(round(L_p1/h)+2:round((L_p1+L_n)/h)+1,1) + N_d);
    rho_v(round((L_p1+L_n)/h)+2:length(X)-1,1) = q*(p(round((L_p1+L_n)/h)+2:length(X)-1,1) - n(round((L_p1+L_n)/h)+2:length(X)-1,1) + N_d1);
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
rho_v(1:round(L_p1/h)+1,1) = q*(p(1:round(L_p1/h)+1,1) - n(1:round(L_p1/h)+1,1) - N_a1);
rho_v(round(L_p1/h)+2:round((L_p1+L_n)/h)+1,1) = q*(p(round(L_p1/h)+2:round((L_p1+L_n)/h)+1,1) - n(round(L_p1/h)+2:round((L_p1+L_n)/h)+1,1) + N_d);
rho_v(round((L_p1+L_n)/h)+2:length(X),1) = q*(p(round((L_p1+L_n)/h)+2:length(X),1) - n(round((L_p1+L_n)/h)+2:length(X),1) + N_d1);
rho_v = rho_v';

% Potential in volts
V = V';
V2 = V;

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

L_n = 1; % Length n region in micrometer

h = 0.01; % Step size in mirometer
X = -L_p1:h:(L_n+L_n1); % X-axis
X1 = X;

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

for k = 1:1000

    V(1,1) = -V_t*log(N_a1/n_i); % Boundary Condition for n-plus region end
    V(length(X),1) = V_t*log(N_d1/n_i); % Boundary Condition for n region end
    
    p = n_i * exp(-V/V_t);
    n = n_i * exp(V/V_t);

    % Charge Density Matrix
    rho_v(1,1) = 0; % Boundary Condition for n-plus region end
    rho_v(2:round(L_p1/h)+1,1) = q*(p(2:round(L_p1/h)+1,1) - n(2:round(L_p1/h)+1,1) - N_a1);
    rho_v(round(L_p1/h)+2:round((L_p1+L_n)/h)+1,1) = q*(p(round(L_p1/h)+2:round((L_p1+L_n)/h)+1,1) - n(round(L_p1/h)+2:round((L_p1+L_n)/h)+1,1) + N_d);
    rho_v(round((L_p1+L_n)/h)+2:length(X)-1,1) = q*(p(round((L_p1+L_n)/h)+2:length(X)-1,1) - n(round((L_p1+L_n)/h)+2:length(X)-1,1) + N_d1);
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

% Potential in volts
V = V';
V1 = V;

L_n = 3; % Length n region in micrometer

h = 0.01; % Step size in mirometer
X = -L_p1:h:(L_n+L_n1); % X-axis
X3 = X;

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

for k = 1:1000

    V(1,1) = -V_t*log(N_a1/n_i); % Boundary Condition for n-plus region end
    V(length(X),1) = V_t*log(N_d1/n_i); % Boundary Condition for n region end
    
    p = n_i * exp(-V/V_t);
    n = n_i * exp(V/V_t);

    % Charge Density Matrix
    rho_v(1,1) = 0; % Boundary Condition for n-plus region end
    rho_v(2:round(L_p1/h)+1,1) = q*(p(2:round(L_p1/h)+1,1) - n(2:round(L_p1/h)+1,1) - N_a1);
    rho_v(round(L_p1/h)+2:round((L_p1+L_n)/h)+1,1) = q*(p(round(L_p1/h)+2:round((L_p1+L_n)/h)+1,1) - n(round(L_p1/h)+2:round((L_p1+L_n)/h)+1,1) + N_d);
    rho_v(round((L_p1+L_n)/h)+2:length(X)-1,1) = q*(p(round((L_p1+L_n)/h)+2:length(X)-1,1) - n(round((L_p1+L_n)/h)+2:length(X)-1,1) + N_d1);
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

% Potential in volts
V = V';
V3 = V;

figure;
plot(X1*10^(-4),V1,'Displayname','L_n = 1 um ');
hold on
plot(X2*10^(-4),V2,'Displayname','L_n = 2 um ');
plot(X3*10^(-4),V3,'Displayname','L_n = 3 um ');
xlabel('x ( cm )');
ylabel('V ( volt )');
title("Potential Profile");
legend;
grid;
hold off;