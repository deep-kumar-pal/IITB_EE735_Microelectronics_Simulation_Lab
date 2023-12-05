% Problem 1 : Abrupt junction %
% =========================== %

a = 3; % second last digit
b = 5; % last digit
L_n = 1.5; % Length of n region in micrometer
L_p = 1.5; % Length of p region in micrometer
epsilon_Si = 11.8; % Dielectric constant of Silicon
epsilon_0 = 8.85 * 10^(-14); % Permittivity of free space in F/cm
n_i = 1.5 * 10^10; % Intrinsic carrier concentration per cm^3
N_a = (1+a) * 10^15; % Acceptor doping concentration per cm^3
N_d = (1+b) * 10^15; % Donor doping concentration per cm^3
T = 300; % Temperature in K
q = 1.6*10^(-19); % Charge of an electron
V_t = 0.026; % Thermal voltage in volts

%%% With Depletion Approximation %%%

V_bi = V_t*log((N_a*N_d)/n_i^2);
W_d = sqrt(2*epsilon_Si*epsilon_0*V_bi*q^(-1)*(N_a^(-1)+N_d^(-1))); % Depletion width in cm
W_n = 10^4*(N_a/(N_a+N_d))*W_d; % Depletion width in n-side in micrometer
W_p = 10^4*(N_d/(N_a+N_d))*W_d; % Depletion width in p-side in micrometer

h_test = 0.0001; % Step size
X_test = -L_p:h_test:L_n; % X-axis

% Volume Charge Density matrix in C/cm^3
rho_v_test = zeros(1,length(X_test));
rho_v_test(1,round((-W_p+L_p)/h_test+1):round((L_p/h_test)+1)) = (-q)*N_a;
rho_v_test(1,round((L_p/h_test)+2):round((W_n+L_p)/h_test+1)) = q*N_d;

% Electric Field Profile
E_test = zeros(1,length(X_test));

for i = 2:length(X_test)
    E_test(1,i) = E_test(1,i-1) + h_test * 0.5 * (rho_v_test(1,i-1) + rho_v_test(1,i));
end

E_test = 10^(-4) * (epsilon_Si * epsilon_0)^(-1) * E_test;

% Voltage Profile
V_test = zeros(1,length(X_test));

for i = 2:length(X_test)
    V_test(1,i) = V_test(1,i-1) + h_test * 0.5 * (E_test(1,i-1) + E_test(1,i));
end

V_test = -(10^(-4))*V_test;

%%% Without Depletion Approximation %%%

h = 0.01; % Step size in mirometer
X = -L_p:h:L_n; % X-axis

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

    V(1,1) = -V_t*log(N_a/n_i); % Boundary Condition for p-region end
    V(length(X),1) = V_t*log(N_d/n_i); % Boundary Condition for n-region end

    p = n_i * exp(-V/V_t);
    n = n_i * exp(V/V_t);
    
    % Charge Density Matrix
    rho_v(1,1) = 0; % Boundary Condition for p-region end
    rho_v(2:round((L_p/h)+1),1) = q*(p(2:round((L_p/h)+1),1) - n(2:round((L_p/h)+1),1) - N_a);
    rho_v(round((L_p/h)+2):length(X)-1,1) = q*(p(round((L_p/h)+2):length(X)-1,1) - n(round((L_p/h)+2):length(X)-1,1) + N_d);
    rho_v(length(X),1) = 0; % Boundary Condition for n-region end
    
    % Newton Raphson using Jacobian
    B = -(h^2*10^(-8)*(epsilon_Si*epsilon_0)^(-1))*rho_v;
    f = A * V - B;
    
    for i = 1:length(X)
        B_1(i,i) = -q*(h^2*10^(-8)*(epsilon_Si*epsilon_0)^(-1))*(V_t)^(-1)*n_i*(-exp(-V(i,1)/V_t)-exp(V(i,1)/V_t));
    end

    J = A - B_1;

    % L-U Decomposition
    L = eye(length(X));
    U = J;
    for i = 2:length(X)
        for j = i:length(X)
            L(j,i-1) = U(j,i-1)/U(i-1,i-1);
            U(j,:) = U(j,:) - (U(j,i-1)/U(i-1,i-1))*U(i-1,:);
        end
    end

    y = zeros(length(X),1);
    for i = 1:length(X)
        a = 0;
        for j = 1:length(X)
            if(i~=j)
                a = a + L(i,j)*y(j,1);
            end
        end
        y(i,1) = f(i,1) - a;
    end
    
    x = zeros(length(X),1);
    for i = length(X):-1:1
        a = 0;
        for j = 1:length(X)
            if(i~=j)
                a = a + U(i,j)*x(j,1);
            end
        end
        x(i,1) = (y(i,1) - a)/U(i,i);
    end
    
    V = V - x;

end

% Charge Density in C/cm^3
p = n_i * exp(-V/V_t);
n = n_i * exp(V/V_t);
rho_v(1:round((L_p/h)+1),1) = q*(p(1:round((L_p/h)+1),1) - n(1:round((L_p/h)+1),1) - N_a);
rho_v(round((L_p/h)+2):length(X),1) = q*(p(round((L_p/h)+2):length(X),1) - n(round((L_p/h)+2):length(X),1) + N_d);
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
plot(X_test* 10^(-4),V_test,'Displayname','with approx.');
hold on;
plot(X* 10^(-4),V-V(1,1),'Displayname','without approx.');
xlabel('x ( cm )');
ylabel('V ( volt )');
title('Voltage Profile');
legend;
grid on;
hold off;

figure;
plot(X_test* 10^(-4),E_test,'Displayname','with approx.');
hold on;
plot(X* 10^(-4),E,'Displayname','without approx.');
xlabel('x ( cm )');
ylabel('E ( V/cm )');
title('Electric Field Profile');
legend;
grid on;
hold off;

figure;
plot(X_test* 10^(-4),rho_v_test/q,'Displayname','with approx.');
hold on
plot(X* 10^(-4),rho_v/q,'Displayname','without approx.');
xlabel('x ( cm )');
ylabel('Charge Concentration ( /cm^3 )');
title('Charge Concentration Profile');
legend;
grid on;
hold off;

fprintf("\nBuilt-in potential (with depletion approximation) = %f V",V_test(1,length(X_test))-V_test(1,1));
fprintf("\nBuilt-in potential (without depletion approximation) = %f V",V(1,length(X))-V(1,1));

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
grid on;
hold off;

% Energy Band Diagram

E_F = zeros(1,length(X));
E_i = E_F - V;
E_g = 1.12; 
E_C = E_i + E_g/2;
E_V = E_i - E_g/2;

figure;
plot(X,E_C,'Displayname','E_C');
hold on
plot(X,E_i,"--",'Displayname','E_{mid}');
plot(X,E_F,"k",'Displayname','E_F');
plot(X,E_V,"r",'Displayname','E_V');
h = gca;
h.XAxis.Visible = 'off';
ylabel("(E - E_F) in eV");
title("Energy Band Diagram");
legend;
grid on;
hold off;