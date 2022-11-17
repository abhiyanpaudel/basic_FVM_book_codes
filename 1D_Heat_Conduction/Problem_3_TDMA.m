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
T(1) = T0;

% Initialization
for i = 2:n
T(i) = 170;
end


nIter = 0;
ErrorMax = 10^-4;
FCmax = 0.1;
st = ['A' 'T'];
fid = fopen('Temperature_dist_Annular_Fin_TDMA.txt','w');
fprintf(fid,'x     0 mm        2.083 mm    6.25 mm    10.417 mm    14.58 mm    18.75 mm    22.917 mm   25.0 mm\n');
 rr = r2-r1;
for i = 1:n
    if x(i) < rr
        k(i) = k2;
    else
        k(i) = k3;
    end
end
for i = 1:n
    C = 2*pi*(r1+xcf(i));
    A(i) = C*t;
    Peri(i) = 4*pi*(r1+x(i));
end
F = 10^6*[A(1);A(2);A(3);A(4);A(5);A(6);A(7);A(n)];
 fprintf(fid,'%s%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f\n',st(1),F);
for i = 2:n-1
SP(i) = h*Peri(i)*(xcf(i+1)-xcf(i));
SU(i) = SP(i)*T_inf;
end

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
%  CALCULATE CELL FACE CONDUCTIVITY BY HARMONIC MEAN.       
    dxe = x(i+1)-x(i);
    dxep = x(i+1)-xcf(i+1);
    dxem = xcf(i+1)-x(i);
    dxw = x(i)-x(i-1);
    dxwp = x(i)-xcf(i);
    dxwm = xcf(i)-x(i-1);
    ksme = dxe/(dxem/k(i)+dxep/k(i+1))*(1-LE)+LE*k(i+1);
    ksmw = dxw/(dxwp/k(i)+dxwp/k(i-1))*(1-LW)+LW*k(i-1);
    AW(i) = ksmw*A(i)/dxw;
    AE(i) = ksme*A(i+1)/dxe;
end

%TriDiagonal Matrix Algorithm (TDMA) or Thomas Algorithm
% a_i*T_i = b_i*T_(i+1) + c_i*T_(i-1) + d_i
% a,b,c,d are input vectors. T is the solution, also a vector.

% Performs Gauss elimination
for i = 2:n-1
    AP(i) = AE(i)+AW(i);
    a = AP(i)+SP(i);
    b = AE(i);
    c = AW(i);
    d = SU(i);
    if i == 2
        P(i) = b/a;
        Q(i) = (d+c*T(1))/a;
    else
        P(i) = b/(a-c*P(i-1));
        Q(i) = (d+c*Q(i-1))/(a-c*P(i-1));
    end
    
end
while FCmax>ErrorMax
    nIter = nIter+1;
    Told = T;
    % Backward substitution
    for i = n-1:-1:2
        if i == n-1
            T(i) = Q(i)/(1-P(i));
        else
            T(i) = P(i)*T(i+1)+Q(i);
        end
    end
    T(n) = T(n-1);
   
    for i = 1:n
        FC(i) = (T(i)-Told(i))/Told(i);
    end
    FCmax = max(FC)
  
end
  fprintf(fid,'%s%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f\n',st(2),T(1),T(2),T(3),T(4),T(5),T(6),T(7),T(n));
Q = -k2*A(1)*(T(2)-T(1))/(x(2)-x(1))               % Heat loss from the fin
plot(x,T,'-*')
xlabel('x')
ylabel('T')
