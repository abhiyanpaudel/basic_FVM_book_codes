clc
clear all
close all

% Parameters
l = 0.02;                              % Length of a fin
t = 0.002;                             % Thickness of a fin
b = 0.2;                               % Breadth of a fin
T_w = 225;                             % Wall temperature;
T_inf = 25;                            % Ambient Temperature
k = 45;                                % Thermal conductivity of fin material
h = 15;                                % Heat transfer coefficient
A = b*t;                               % Cross sectional area of a fin
P = 2*b;                               % Perimeter of a fin
dx = 0.004;                            % Spatial step
n = 7;                                 % No. of nodes

% Cell face coordinate 
xcf(1) = 0;
xcf(n) = l;
for i = 2:n-1
    if i == 2
        xcf(i) = xcf(1);
    else
        xcf(i) = xcf(2)+(i-2)*dx;
    end
end

% Node coordinate
x(1) = xcf(1);
for i = 2:n-1
    x(i) = 0.5*(xcf(i)+xcf(i+1));
end
x(n) = xcf(n);
% Boundary condition
T(1) = T_w;

% Thermal conductivity of fin material
for i = 1:n
    k(i) = 45;      
end

% Initialization

for i = 2:n
    if i == 2 ||  i == 7
        T(i) = T(i-1)-2;
    else
        T(i) = T(i-1)-4;
    end
    
end

Spi = h*P*dx;
Sui = Spi*T_inf;
nIter = 0;
ErrorMax = 10^-4;
FCmax = 0.1;
fid = fopen('Temperature_dist_Rectangular_Fin.txt','w');
fprintf(fid,'l       FCMAX      0 cm        0.2 cm     0.6 cm      1.0 cm      1.4 cm      1.8 cm      2 cm\n');
fprintf(fid,'%d            %12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f\n',nIter,T(1),T(2),T(3),T(4),T(5),T(6),T(n));
% coefficients

for i = 2:n-1
    if i == 2 
        LW = 1;
    else
        LW = 0;
    end
    
    if i == n-1
        LE = 1;
    else
        LE = 0;
    end
        
    dxe = x(i+1)-x(i);
    dxep = x(i+1)-xcf(i+1);
    dxem = xcf(i+1)-x(i);
    dxw = x(i)-x(i-1);
    dxwp = x(i)-xcf(i);
    dxwm = xcf(i)-xcf(i-1);
    ksme = dxe/(dxem/k(i)+dxep/k(i+1))*(1-LE)+LE*k(i+1);
    ksmw = dxw/(dxwp/k(i)+dxwp/k(i-1))*(1-LW)+LW*k(i-1);
    AW(i) = ksmw*A/dxw;
    AE(i) = ksme*A/dxe;
end

while FCmax>ErrorMax
    nIter = nIter+1;
    Told = T;
    for i = 2:n-1
        AP(i) = AE(i)+AW(i);
        if i == n-1
            T(i) = (AW(i)*T(i-1)+Sui)/(Spi+AP(i)-AE(i));             % Boundary condition, q = 0 at n = 7 such that T6 = T7;
        else
            T(i) = (AE(i)*T(i+1)+AW(i)*T(i-1)+Sui)/(Spi+AP(i));
        end
    end
    T(n) = T(n-1);
    for i = 2:n-1
        FC(i) = abs(T(i)-Told(i))/Told(i);
    end
    FCmax = max(FC)
    fprintf(fid,'%d%12.7f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f\n',nIter,FCmax,T(1),T(2),T(3),T(4),T(5),T(6),T(n));
end
Qloss = AW(2)*(T(1)-T(2))                                            % fin heat loss
% Exact Solution
m = sqrt(h*P/(k(1)*A));
for i = 2:n
    T(i) = ((T_w-T_inf)*cosh(m*(l-x(i)))/cosh(m*l))+T_inf;
end
st = 'Exact';
fprintf(fid,'%s         %12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f\n',st,T(1),T(2),T(3),T(4),T(5),T(6),T(n));

fprintf('Total number of iterations = %d\n',nIter)


