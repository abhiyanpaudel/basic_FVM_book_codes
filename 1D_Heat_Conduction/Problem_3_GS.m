clc
clear all
close all

% Parameters
k2 = 200;                                    % Thermal conductivity of inner material(W/m-K)
k3 = 40;                                     % Thermal conductivity of outer material(W/m-K)
r1 = 0.0125;                                 % Radius of pipe(m) 
r2 = 0.025;                                  % Radius of inner material(m)
r3 = 0.0375;                                 % Radius of outer material(m)   
t = 0.001;                                   % Fin thickness(m)             
h = 20;                                      % Heat transfer coefficient(W/m^2-K)
T0 = 200;                                    % Fin base temperature(C)
T_inf = 25;                                  % Ambient temperature(C)
n = 8;                                       % Total number of nodes;
l = r3-r1;                                   % Distance from fin base to fin tip(m)
dx = (r3-r1)/(n-2);                          % Spatial step


% Discretization
x(1) = 0;
x(n) = l;
for i = 2:n-1
    if i == 2
        x(i) = x(1)+dx/2;
    else
        x(i) = x(2)+(i-2)*dx;
    end
end

% Boundary condition
T(1) = T0;                                   

for i = 2:n
    T(i) = 170;
end

nIter = 0;
ErrorMax = 10^-4;
FCmax = 0.1;
fid = fopen('Temperature_dist_Annular_Fin_GS.txt','w');
fprintf(fid,'l       FCMAX      0 mm        2.083 mm    6.25 mm    10.417 mm    14.58 mm    18.75 mm    22.917 mm   25.0 mm\n');
fprintf(fid,'%d            %12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f\n',nIter,T(1),T(2),T(3),T(4),T(5),T(6),T(7),T(n));

for i = 1:n
     C = 2*pi*(r1+x(i));
     A(i) = C*t;
     P(i) = 2*C;
end
while FCmax>ErrorMax
    nIter = nIter+1;
    Told = T;
    for i = 2:n-1
        Spi = h*P(i)*dx;
        Sui = Spi*T_inf;
        if (x(i)+r1) < r2
        AE = (k2*A(i))/(x(i+1)-x(i));
        AW = (k2*A(i))/(x(i)-x(i-1));
        else
        AE = (k3*A(i))/(x(i+1)-x(i));
        AW = (k3*A(i))/(x(i)-x(i-1)); 
        end
        AP = AE+AW;
        if i == n-1
            T(i) = (AW*T(i-1)+Sui)/(Spi+AP-AE);             % Boundary condition, q = 0 at n = 8 such that T7 = T8;
        else
            T(i) = (AE*T(i+1)+AW*T(i-1)+Sui)/(Spi+AP);
        end
    end
    T(n) = T(n-1);
    for i = 1:n
        FC(i) = (T(i)-Told(i))/Told(i);
    end
    FCmax = max(FC);
    fprintf(fid,'%d%12.7f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f\n',nIter,FCmax,T(1),T(2),T(3),T(4),T(5),T(6),T(7),T(n));
end


plot(x,T,'*-')
