function x=backsub(A,b)
%Lb=b_prime;
%Ux=b;

% This function performs back substitution on an upper triangular matrix that has
% been modified by concatenating the RHS of the system.  
% Note that B is assumed to be upper triangular at this point.
Aback=cat(2,A,b);
n=size(Aback,1);                   %number of unknowns in the system
x=zeros(n,1);          %space in which to store our solution vector
x(n)=Aback(n,n+1)/Aback(n,n);          %finalized solution for last variable, resulting from upper triangular conversion

for ir1=n-1:-1:1
    x(ir1)=Aback(ir1,n+1);       %assume we're only dealing with a single right-hand side here.
    fact=Aback(ir1,ir1);         %diagonal element to be divided through doing subs for the ir2 row
    for ic=ir1+1:n
        x(ir1)=x(ir1)-Aback(ir1,ic)*x(ic);
    end %for
    x(ir1)=x(ir1)/fact;      %divide once at the end to minimize number of ops
end %for

end %function