clc
clear all
close all

% Parameters 
 k = 0.25;                                    % Thermal conductivity(W/m-k)
 C = 2000;                                    % Specific heat(J/kg-k)
 rho = 1300;                                  % Density of plastic sheet(kg/m^3)
 A = 1;                                       % Area(m^2)
 n = 7;                                       % Total number of nodes;
 l = 0.01;                                    % Total thickness of two plastics   
 dx = 0.002;                                  % Spatial step
 dt = 20;                                     % Time step
 
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
 
 % Initialization 
 
 for i = 2:n-1
     T(i) = 30;
 end
 
 % Boundary condition 
 Tbound = 250;
 T(1) = Tbound;
 T(n) = Tbound;
 
 Tmid = 30;
 cf = rho*A*dx*C/dt;
 t = 0;
 fid = fopen('Temperature_dist_explicit.txt','w');
%fprintf(fid,'% s%8s%8s%8s%8s%8s%8s%8s\n','Time','0mm','1mm','3mm','5mm','7mm','9mm','10mm');
fprintf(fid,'Time    0 mm      1 mm       3 mm      5 mm      7 mm     9 mm      10 mm\n');
fprintf(fid,'% d%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f\n',t,T(1),...
    T(2),T(3),T(4),T(5),T(6),T(n));
while Tmid<140
    t = t+dt;
    Tnew(1) = Tbound;
    for i = 2:n-1 
        AE = (k*A)/(x(i+1)-x(i));
        AW = (k*A)/(x(i)-x(i-1));
        Tnew(i) = (AE*T(i+1)+AW*T(i-1)+(cf-(AE+AW))*T(i))/cf;   % By equation 1.31
    end 
    Tnew(n) = Tbound;
    T = Tnew;
    Tmid = T(4);
    if t <99
    fprintf(fid,'% d%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f\n',t,T(1),T(2),...
        T(3),T(4),T(5),T(6),T(n));
    else
      fprintf(fid,'%d%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f\n',t,T(1),T(2),...
        T(3),T(4),T(5),T(6),T(n));
    end
end
 fprintf('Time required to press two sheets = %5.2f\n',t)
 
plot(x,T,'-ro','LineWidth',3,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','g',...
                       'MarkerSize',10)
xlabel('x')
ylabel('T')
 title('Temperature Distribution over the two plastic sheets')
 
 
 
 
 