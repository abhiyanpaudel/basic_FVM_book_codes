clc
clear all
close all

% Problem Data

rho = 1;                                      % Density of a fluid (kg/m^3)
L = 2;                                        % Nozzle length (m)
xa = 0;                                       % Position of first node A(m)
xb = L;                                       % Position of last node E(m)
m = 4;                                        % Number of control volumes 
A_A = 0.5;                                    % Cross Sectional area at the inlet(m^2)
A_E = 0.1;                                    % Cross Sectional area at the outlet(m^2)
P0 = 10;                                      % Stagnation Pressure at inlet(Pa)
PE = 0;                                       % Static pressure at exit(Pa)
m_dot = 1;                                    % Mass flow rate (kg/s)
deltax = abs(xb-xa)/m;                        % Grid size
xnodes = xa:deltax:xb;                        % Position of main grids(m)
xcenters = (xa+deltax/2):deltax:(xb-deltax/2);% Postion of volume centers(m)


% Compute area at nodes and cell centers 
% At nodes
for i = 1:length(xnodes)
    nodeA(i) = A_E+2*0.1*(2-xnodes(i));
end
% At cell centres
for i = 1:length(xcenters)
centerA(i) = A_E+2*0.1*(2-xcenters(i));
end

% Velocity initialization
u = zeros(1,length(xcenters));

% Compute guessed velocity using constant mass flow rate
for i = 1:length(xcenters)
    u(i) = m_dot/(rho*centerA(i));
end
P = zeros(1,length(xnodes));
% Guessed values of pressure at nodes
for i = 1:length(xnodes)
    P(i) = P0-(i-1)*2.5;
end

FCmax = 0.1;
Errormax = 10^-7;
nIter = 0;

while FCmax>Errormax
    nIter = nIter+1
    uold = u;
    Pold = P;
%Discretized u-momentum equation 
for i = 1:length(xcenters)
     if i == 1
          uA = u(i)*centerA(i)/nodeA(i);
          Fw = rho*uA*nodeA(i);
          Fe = rho*(u(i)+u(i+1))*0.5*nodeA(i+1);
          aW(i) = 0;
          aE(i) = 0;
          aP(i) = Fe+Fw*0.5*(centerA(i)/nodeA(i))^2;
          Su(i)= (P0-P(i+1))*centerA(i)+Fw*(centerA(i)/nodeA(i))*u(i);
          d(i) = centerA(i)/aP(i);
     elseif i == length(xcenters)
          Fw = rho*(u(i-1)+u(i))*0.5*nodeA(i);
          Fe = rho*u(i)*centerA(i);
          aW(i) = Fw;
          aE(i) = 0;
          aP(i) = aW(i)+aE(i)+(Fe-Fw);
          Su(i)= (P(i)-P(i+1))*centerA(i);
          d(i) = centerA(i)/aP(i);
     else
          Fw = rho*(u(i-1)+u(i))*0.5*nodeA(i);
          Fe = rho*(u(i)+u(i+1))*0.5*nodeA(i+1);
          aW(i) = Fw;
          aE(i) = 0;
          aP(i) = aW(i)+aE(i)+(Fe-Fw);
          Su(i)= (P(i)-P(i+1))*centerA(i);
          d(i) = centerA(i)/aP(i);
     end
 end


 %Forward substitution 
% p_i*T_i = q_i*T_(i+1) + r_i*T_(i-1) + s_i
% p,q,r,s are input vectors. T is the solution, also a vector.

for i = 1:length(xcenters)
    p = aP(i);
    q = aE(i);
    r = aW(i);
    s = Su(i);
    if i == 1
        u(i) = s/p;
    else
        u(i) = (r*u(i-1)+s)/p;
    end
    
end

P_prime = zeros(1,length(xnodes));
% Pressure Correction Equation

for i = 2:length(xnodes)-1
    AW(i) = rho*d(i-1)*centerA(i-1);
    AE(i) = rho*d(i)*centerA(i);
    FW(i) = rho*u(i-1)*centerA(i-1);
    FE(i) = rho*u(i)*centerA(i);
    AP(i) = AW(i)+AE(i);
    b_prime(i) = FW(i)-FE(i);
end

 %TriDiagonal Matrix Algorithm (TDMA) or Thomas Algorithm
% p_i*T_i = q_i*T_(i+1) + r_i*T_(i-1) + s_i
% p,q,r,s are input vectors. Phi is the solution, also a vector.

% Performs Gauss elimination
for i = 2:length(xnodes)-1
    p = AP(i);
    q = AE(i);
    r = AW(i);
    s = b_prime(i);
    if i == 2
        A(i) = q/p;
        B(i) = (s+r*P_prime(i-1))/p;
    else
        A(i) = q/(p-r*A(i-1));
        B(i) = (s+r*B(i-1))/(p-r*A(i-1));
    end
    
end
% Performs backward substitution
for i = length(xnodes)-1:-1:2
  P_prime(i) = A(i)*P_prime(i+1)+B(i);
end

%Corrected nodal pressures using pressure corrections
for i = 2:length(xnodes)-1
    P(i) = P(i)+P_prime(i);
end


%Corrected velocities
for i = 1:length(xcenters)
    u(i) = u(i)+d(i)*(P_prime(i)-P_prime(i+1));
end


% Corrected nodal pressure at A
P(1) = P0-0.5*rho*u(1)^2*(centerA(1)/nodeA(1))^2;

% Continuity check
for i = 1:length(xcenters)
    mass_rate(i) = rho*u(i)*centerA(i);
end
% Error monitor
    for i = 1:length(xcenters)
        FC(i) = abs((u(i)-uold(i))/uold(i));
    end
    FCmax = max(FC)
end


plot(xnodes,0,'-*r')
hold on 
plot(xcenters,0,'*b')
