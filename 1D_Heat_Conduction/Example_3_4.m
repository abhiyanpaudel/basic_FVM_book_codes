                              % flow velocity(m/s)
 
clc
clear all
close all

% Parameters 

 L = 1;                                       % Length of a fin(m)  
 rho = 1;                                     % Density(kg/m^3)
 gamma = 0.1;                                 % Diffusion Coefficient
 m = 7;                                      % Total number of nodes
 dx = L/(m-2);                                % Grid size
 Phi0 = 1;                                    % Phi at x = 0
 PhiL = 0;                                    % Phi at x = L
 u = 0.2;                                     % Flow velocity
 F = rho*u;                                   % Convective mass flux per unit area 
 D = gamma/dx;                                % Diffusion conductance
 Pe = F/D;                                     % Peclet Number 
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
         AWW(i) = 0;
         AW(i) = 0;
         AE(i) = D+D/3-(3/8)*F;
         SU(i) = ((8/3)*D+(2/8)*F+F)*Phi(i-1);
         SP(i) = -((8/3)*D+(2/8)*F+F);
         AP(i) = AWW(i)+AW(i)+AE(i)-SP(i);
     elseif i == m-1
         AWW(i) = -(1/8)*F;
         AW(i) = D+(1/3)*D+(6/8)*F;
         AE(i) = 0;
         SU(i) = ((8/3)*D-F)*Phi(i+1);
         SP(i) = -((8/3)*D-F);
         AP(i) = AWW(i)+AW(i)+AE(i)-SP(i);
     elseif i == 3
         AWW(i) = 0;
         AW(i) = D+(7/8)*F+(1/8)*F;
         AE(i) = D-(3/8)*F;
         SU(i) = -F*Phi(i-2)/4;
         SP(i) = F/4;
         AP(i) = AWW(i)+AW(i)+AE(i)-SP(i);
     else
         AWW(i) = -F/8;
         AW(i) = D+F*6/8+F*1/8;
         AE(i) = D-F*3/8;
         SU(i) = 0;
         SP(i) = 0;
         AP(i) = AWW(i)+AW(i)+AE(i)-SP(i);
     end
     
 end
 

FCmax = 0.1;
Errormax = 10^-4;
nIter = 0;
% Backward substitution
while FCmax>Errormax
    nIter = nIter+1
    phi_old = Phi;
    for i = 2:m-1
        if i == 2
            Phi(i) = (AW(i)*Phi(i-1)+AE(i)*Phi(i+1)+SU(i))/AP(i);
        else
            Phi(i) = (AW(i)*Phi(i-1)+AE(i)*Phi(i+1)+AWW(i)*Phi(i-2)+SU(i))/AP(i);
        end
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
xlabel('Distance x(m)')
 ylabel('\phi')
legend('Finite Volume Solution','Analytical Solution','Location','NorthEast')
legend boxoff
text(0.8,0.7,'\itu = \rm0.2m/s')


 