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
%      Equation (2.18)
     if i == 2
         AW(i) = 0;
         AE(i) = k*A/dx;
         SU(i) = 2*k*A/dx*T(i-1);
         SP(i) = -2*k*A/dx;
         AP(i) = AE(i)+AW(i)-SP(i);
%      Equation (2.19)
     elseif i == n-1
         AW(i) = k*A/dx;
         AE(i) = 0;
         SU(i) = 2*k*A/dx*T(i+1);
         SP(i) = -2*k*A/dx;
         AP(i) = AE(i)+AW(i)-SP(i);
%      Equation (2.17)
      else
         AW(i) = k*A/dx;
         AE(i) = k*A/dx;
         SU(i) = 0;
         SP(i) = 0;
         AP(i) = AE(i)+AW(i)-SP(i);
     end
 end
 
 %TriDiagonal Matrix Algorithm (TDMA) or Thomas Algorithm
% a_i*T_i = b_i*T_(i+1) + c_i*T_(i-1) + d_i
% a,b,c,d are input vectors. T is the solution, also a vector.

% Performs Gauss elimination
for i = 2:n-1
    a = AP(i);
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

FCmax = 0.1;
Errormax = 10^-4;
nIter = 0;
% Backward substitution
while FCmax>Errormax
    nIter = nIter+1;
    Told = T;
    for i = n-1:-1:2
            T(i) = P(i)*T(i+1)+Q(i);
    end
    for i = 2:n-1
        FC(i) = (T(i)-Told(i))/Told(i);
    end
    FCmax = max(FC);
end

% Print AW,AE,AP,SU,SP
fprintf('Node    AW      AE        SU        SP        AP\n')
for i = 2:n-1
    fprintf('%d%10.2f%10.2f%10.2f%10.2f%10.2f\n',i,AW(i),AE(i),SU(i),SP(i),AP(i));
end

% Exact Solution
% T = ((TB-TA)/L)*x+TA
T_exact(1) = TA;
T_exact(n) = TB;
for i = 2:n-1
    T_exact(i) = ((TB-TA)/l)*x(i)+TA;
end
 plot(x,T,'sb')
 hold on
 plot(x,T_exact,'k-','LineWidth',2)
 axis([0 0.5 0 500])
 xlabel('Distance x(m)')
 ylabel('Temperature ( ^{o}C)')
legend('Finite Volume Solution','Analytical Solution','Location','East')
legend boxoff

 