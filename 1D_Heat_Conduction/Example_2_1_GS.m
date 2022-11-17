clc
clear all
close all

% Parameters 
 k = 1000;                                    % Thermal conductivity(W/m-k)
 A = 10*10^-3;                                % Area(m^2)
 n = 7;                                       % Total number of nodes;
 l = 0.5;                                     % Length of a rod(m)   
 dx = 0.1;                                    % Spatial step
 TA = 100;                                    % Temperature at A in degree celsius
 TB = 500;                                    % Temperature at B in degree celsius
 
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
 
%  Boundary Condition
T(1) = TA;
T(n) = TB;

% Initialization
for i = 2:n-1
    T(i) = 150;
end
 
 for i = 2:n-1
     if i == 2
         AW(i) = 0;
         AE(i) = k*A/dx;
         SU(i) = 2*k*A/dx*T(i-1);
         SP(i) = -2*k*A/dx;
         AP(i) = AE(i)+AW(i)-SP(i);
     elseif i == n-1
         AW(i) = k*A/dx;
         AE(i) = 0;
         SU(i) = 2*k*A/dx*T(i+1);
         SP(i) = -2*k*A/dx;
         AP(i) = AE(i)+AW(i)-SP(i);
      else
         AW(i) = k*A/dx;
         AE(i) = k*A/dx;
         SU(i) = 0;
         SP(i) = 0;
         AP(i) = AE(i)+AW(i)-SP(i);
     end
 end
 nIter = 0;
ErrorMax = 10^-4;
FCmax = 0.1;
 
while FCmax>ErrorMax
    nIter = nIter+1
    Told = T;
     for i = 2:n-1
            T(i) = (AE(i)*T(i+1)+AW(i)*T(i-1)+SU(i))/AP(i);
    end
 
    
    for i = 2:n-1
        FC(i) = (T(i)-Told(i))/Told(i);
    end
    FCmax = max(FC)

end

 plot(x,T,'*-')