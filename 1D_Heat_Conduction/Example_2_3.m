clc
clear all
close all

% Parameters 

 L = 1;                                       % Length of a fin(m)   
 dx = 0.2;                                    % Grid size
 m = 2+L/dx;                                  % Total number of nodes;
 TB = 100;                                    % Base Temperature in degree celsius
 Tinf = 20;                                   % Ambient Temperature in degree celsius
 n = 5;                                       % n^2=hp/kA (m^-2)
                                      
 
 % Discretization 
 x(1) = 0;
 for i = 2:m-1
     if i == 2
         x(i) = x(1)+dx/2;
     else
         x(i) = x(2)+(i-2)*dx;
     end
 end
x(m) = L;

%  Boundary Condition
T(1) = TB;


% Initialization
for i = 2:m
    T(i) = 100;
end
 
 for i = 2:m-1
  % Equation 2.48
     if i == 2
         AW(i) = 0;
         AE(i) = 1/dx;
         SU(i) = n^2*dx*Tinf+T(i-1)*2/dx;
         SP(i) = -n^2*dx-2/dx;
         AP(i) = AE(i)+AW(i)-SP(i);
  % Equation 2.53
     elseif i == m-1
         AW(i) = 1/dx;
         AE(i) = 0;
         SU(i) = n^2*dx*Tinf;
         SP(i) = -n^2*dx;
         AP(i) = AE(i)+AW(i)-SP(i);
  % Equation 2.44       
      else
         AW(i) = 1/dx;
         AE(i) = 1/dx;
         SU(i) = n^2*dx*Tinf;
         SP(i) = -n^2*dx;
         AP(i) = AE(i)+AW(i)-SP(i);
     end
 end
 
 %TriDiagonal Matrix Algorithm (TDMA) or Thomas Algorithm
% a_i*T_i = b_i*T_(i+1) + c_i*T_(i-1) + d_i
% a,b,c,d are input vectors. T is the solution, also a vector.

% Performs Gauss elimination
for i = 2:m-1
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
    nIter = nIter+1
    Told = T;
    for i = m-1:-1:2
            T(i) = P(i)*T(i+1)+Q(i);
    end
    for i = 2:m-1
        FC(i) = (T(i)-Told(i))/Told(i);
    end
    FCmax = max(FC)
end
T(m) = T(m-1);
% Print AW,AE,AP,SU,SP
fprintf('Node    AW      AE        SU        SP        AP\n')
for i = 2:m-1
    fprintf('%d%10.2f%10.2f%10.2f%10.2f%10.2f\n',i,AW(i),AE(i),SU(i),SP(i),AP(i));
end


T_exact(1) = TB;
% Exact Solution
% T = Tinf+(TB-Tinf)*(cosh(n(L-x))/cosh(nL))
for i = 2:m-1
    T1 = TB-Tinf;
    L1 = n*(L-x(i));
    T_exact(i) = Tinf+T1*(cosh(L1)/cosh(n*L));
end
T_exact(m) = T_exact(m-1);     % since q = 0

% Comparison between finite volume solution and the exact solution
for i = 2:m-1
    Error_perc(i) = ((T_exact(i)-T(i))/T_exact(i))*100;
end

 plot(x,T,'sb')
 hold on
 plot(x,T_exact,'k-','LineWidth',3)
 xlabel('Distance x(m)')
 ylabel('Temperature ( ^{o}C)')
 legend('Finite Volume Solution','Analytical Solution','Location','NorthEast')
 legend boxoff
  %axis([0 1 20 100])


 