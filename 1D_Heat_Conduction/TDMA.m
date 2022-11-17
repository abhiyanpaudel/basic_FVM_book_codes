function X = TDMA(p,q,r)
%TriDiagonal Matrix Algorithm (TDMA) or Thomas Algorithm
% X_i = p_i*X_(i+1) + q_i*X_(i-1) + r_i
% a,b,c are input vectors. X is the solution, also a vector. 

N = length(p)
X = zeros(N,1);
% Performs Gauss Elimination
A(2) = p(2);
B(2) = r(2);
for i = 3:N+1
dr = (1-q(i)*A(i-1));
A(i) = p(i)/dr;
B(i) = (q(i)*B(i-1)+r(i))/dr;
end
% Backward substitution,
A(N+1) = 0;
for i = N+1:-1:2
    if i == N+1 
        X(i) = B(i)/(1-A(i));
    else
      X(i) = A(i)*X(i+1)+B(i) ;
    end 
end
