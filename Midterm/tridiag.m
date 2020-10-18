%%Problem #1a
%%Write a MATLAB or Python function tridiag() that solves a tridiagonal system of equations
%%using the Thomas algorithm. Verify your solution by applying it to the iterative test problem for
%%HW2 in the EP501 Assignments repository

function [Atom] = tridiag(A,b)
%forward elimination
nref=length(b);                %system size for reference problem
Atom=cat(2,A,b); 
for ir1=2:nref
    em=Atom(ir1,1)./Atom(ir1-1,2);
    Atom(ir1,1)=em;
    Atom(ir1,2)=Atom(ir1,2)-em.*Atom(ir1-1,3);
    b(ir1)=b(ir1)-Atom(ir1,1).*b(ir1-1);
end %for

%% Illustrate back substitution on B using provided Matlab function
x(nref)=b(nref)/Atom(nref,2);
for ir1=nref-1:1:-1
    x(ir1)=(b(ir1)-Atom(ir1,ic)).*x(ic+1)./Atom(ir1,2);
end %for
disp('Elimination/back sub solution:  ');
disp(x);

% xsoln=backsub(Atom);
% disp('Elimination/back sub solution:  ');
% disp(xsoln);
% disp('Matlab,GNU/Octave built-in solution:  ');
% disp(a\b);
end %function 