%%Problem #2
% Solve a diagonally dominant system using Jacobi iteration with successive
% over relaxation

nit=10;
Ait=diag(-1*ones(nit-1,1),-1)+diag(-1*ones(nit-1,1),1)+diag(4*ones(nit,1),0);    %this must be diagonally dominant or else the method won't converge
%Ait=randn(nit,nit);    %see if code can detect non-diagonal dominance and exit gracefully...
x0=randn(nit,1);
bit=ones(nit,1);
tol=1e-9;
disp('Verbose Jacobi iterations:  ')
[xit,iterations]=sorFunc(x0,Ait,bit,tol,true);

disp('Solution with Jacobi iterations:  ')
disp(xit);
disp('Number of iterations required and tolerance:  ')
disp(iterations);
disp(tol);
disp('Matlab built-in solution:  ')
disp(Ait\bit);