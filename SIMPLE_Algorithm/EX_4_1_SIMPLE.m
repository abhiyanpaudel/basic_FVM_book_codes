% This program calculates pressure corrections at main nodes 
% and use them to find velocity at control volume faces
 
clc
clear all
close all

% Problem Data

rho = 1;                                      % Density of a fluid (kg/m^3)
L = 1;                                        % Nozzle length (m)
xa = 0;                                       % Position of first node A(m)
xb = L;                                       % Position of last node E(m)
m =  4;                                       % Number of control volumes 
A = 0.5;                                      % Cross Sectional area at the inlet(m^2)
d = 1;                                        % Ae/aP
U = 10;                                       % Inlet velocity (m/s)
P_D = 0;                                      % Outlet pressure (Pa) 
deltax = abs(xb-xa)/m;                        % Grid size
xfaces = xa:deltax:(xb-deltax/2);             % Position of cell faces(m)
xnodes = (xa+deltax/2):deltax:xb;             % Postion of main nodes(m)

% Velocity initialization
u = zeros(1,length(xfaces));

% Boundary condition 
u(1) = U;                                            

% Pressure initialization
P = zeros(1,length(xnodes));

% Boundary condition 
P(length(xnodes)) = P_D;  

% Compute guessed velocity using constant mass flow rate
for i = 2:length(xfaces)
    u(i) = 8+(i-1)*2;
end

P_prime = zeros(1,length(xnodes));
% Pressure Correction Equation

for i = 1:length(xnodes)-1
   if i == 1
       aW(i) = 0;
       aE(i) = rho*d*A;
       aP(i) = aW(i)+aE(i);
       b_prime(i) = rho*A*(u(i)-u(i+1));
   else
       aW(i) = rho*d*A;
       aE(i) = rho*d*A;
       aP(i) = aW(i)+aE(i);
       b_prime(i) = rho*A*(u(i)-u(i+1));
   end
end

 %TriDiagonal Matrix Algorithm (TDMA) or Thomas Algorithm
% p_i*T_i = q_i*T_(i+1) + r_i*T_(i-1) + s_i
% p,q,r,s are input vectors. Phi is the solution, also a vector.

% Performs Gauss elimination
for i = 1:length(xnodes)-1
    p = aP(i);
    q = aE(i);
    r = aW(i);
    s = b_prime(i);
    if i == 1
        A(i) = q/p;
        B(i) = s/p;
    else
        A(i) = q/(p-r*A(i-1));
        B(i) = (s+r*B(i-1))/(p-r*A(i-1));
    end
    
end
% Performs backward substitution
for i = length(xnodes)-1:-1:1
  P_prime(i) = A(i)*P_prime(i+1)+B(i);
end

%Corrected nodal pressures using pressure corrections
for i = 2:length(xnodes)-1
    P(i) = P(i)+P_prime(i);
end


for i = 2:length(xfaces)
    u(i) = u(i)+d*(P_prime(i-1)-P_prime(i));
end

plot(xnodes,0,'-*r')
hold on 
plot(xfaces,0,'*b')
