                              % flow velocity(m/s)
 
clc
clear all
close all

% Parameters 

 L = 1;                                       % Length of a fin(m)  
 rho = 1;                                     % Density(kg/m^3)
 gamma = 0.1;                                 % Diffusion Coefficient
 m = 7;                                       % Total number of nodes
 dx = L/(m-2);                                % Grid size
 Phi0 = 1;                                    % Phi at x = 0
 PhiL = 0;                                    % Phi at x = L
 u = 2.5;                                     % Flow velocity
 F = rho*u;                                   % Convective mass flux per unit area 
 D = gamma/dx;                                % Diffusion conductance
 
% In this problem area at all faces is taken constant i.e. Ae = Aw =A
% F and D are also taken constant for all faces..
 
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
Phi(1) = Phi0;
Phi(m) = PhiL;


% Initialization
for i = 2:m-1
    Phi(i) = 0.3;
end
 
 for i = 2:m-1
     if i == 2
         AW(i) = 0;
         AE(i) = D+max(0,-F);
         SU(i) = (2*D+F)*Phi(i-1);
         SP(i) = -(2*D+F);
         AP(i) = AE(i)+AW(i)-SP(i);
     elseif i == m-1
         AW(i) = D+max(F,0);
         AE(i) = 0;
         SU(i) = 2*D*Phi(i+1);
         SP(i) = -2*D;
         AP(i) = AE(i)+AW(i)-SP(i);
      else
         AW(i) = D+max(F,0);
         AE(i) = D+max(0,-F);
         SU(i) = 0;
         SP(i) = 0;
         AP(i) = AE(i)+AW(i)-SP(i);
     end
 end
 
 %TriDiagonal Matrix Algorithm (TDMA) or Thomas Algorithm
% a_i*T_i = b_i*T_(i+1) + c_i*T_(i-1) + d_i
% a,b,c,d are input vectors. Phi is the solution, also a vector.

% Performs Gauss elimination
for i = 2:m-1
    a = AP(i);
    b = AE(i);
    c = AW(i);
    d = SU(i);
    if i == 2
        P(i) = b/a;
        Q(i) = (d+c*Phi(1))/a;
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
    phi_old = Phi;
    for i = m-1:-1:2
            Phi(i) = P(i)*Phi(i+1)+Q(i);
    end
    for i = 2:m-1
        FC(i) = (Phi(i)-phi_old(i))/phi_old(i);
    end
    FCmax = max(FC)
end

% % Exact Solution
% % Phi = Tinf+(TB-Tinf)*(cosh(n(L-x))/cosh(nL))
Phi_exact(1) = Phi0;
for i = 2:m-1
    Phi1 = PhiL-Phi0;
    num = exp(F*x(i)/gamma)-1;
    den = exp(F*L/gamma)-1;
    Phi_exact(i) = Phi1*num/den+Phi0;
end
Phi_exact(m) = PhiL;

% Comparison between finite volume solution and the exact solution
for i = 2:m-1
    Error_perc(i) = ((Phi_exact(i)-Phi(i))/Phi_exact(i))*100;
end

 plot(x,Phi,'s-b')
 hold on
 plot(x,Phi_exact,'k-','LineWidth',3)
 axis([0 1 0 1.2])
 xlabel('Distance x(m)')
 ylabel('\phi')
legend('Finite Volume Solution','Analytical Solution','Location','NorthEast')
legend boxoff
text(0.1,0.4,'\itu = \rm2.5m/s')



 