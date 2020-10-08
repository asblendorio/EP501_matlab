function x=backsubwithLinearsolver(L,U,b,b2,b3)
%Solve the system of Linear Equations from test data 
%L_inverse*b=b_prime;
%U*x=b;


Linv=zeros(n,n);
n=length(L);
I=eye(n);
for j = 1:n
    for i = 1:n
        Linv(j,i) = (I(j,i) - L(j,1:(j-1))*Linv(1:(j-1),i))/L(j,j);    % Lower triangular matrix inversion (L^(-1))
    end %for
end %for
% This function performs back substitution on an upper triangular matrix that has
% been modified by concatenating the RHS of the system.  
% Note that B is assumed to be upper triangular at this point.

n=size(A,1);                   %number of unknowns in the system
x=zeros(n,1);                  %space in which to store our solution vector
x(n)=A(n,n+1)/A(n,n);          %finalized solution for last variable, resulting from upper triangular conversion

for ir1=n-1:-1:1
    x(ir1)=A(ir1,n+1);       %assume we're only dealing with a single right-hand side here.
    fact=A(ir1,ir1);         %diagonal element to be divided through doing subs for the ir2 row
    for ic=ir1+1:n
        x(ir1)=x(ir1)-A(ir1,ic)*x(ic);
    end %for
    x(ir1)=x(ir1)/fact;      %divide once at the end to minimize number of ops
end %for

end %function