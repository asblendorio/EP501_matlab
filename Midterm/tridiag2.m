%%Problem #1a
%%Write a MATLAB or Python function tridiag() that solves a tridiagonal system of equations
%%using the Thomas algorithm. Verify your solution by applying it to the iterative test problem for
%%HW2 in the EP501 Assignments repository

function [Atom] = tridiag2(A,b)
%forward elimination
nref=length(b);                %system size for reference problem
Atom=cat(2,A,b); 
Aprime=zeros(nref,3);
for i=2:nref-1
    Aprime(nref,1)=Atom(nref,nref-1);
    Aprime(nref,2)=Atom(nref,nref);
end %for    
Aprime(1,2)=Atom(1,2);
Aprime(1,3)=Atom(1,3);

for ir1=2:nref
    em=Atom(ir1,1)./Atom(ir1-1,2);
    Atom(ir1,1)=em;
    Atom(ir1,2)=Atom(ir1,2)-em.*Atom(ir1-1,3);
    b(ir1)=b(ir1)-Atom(ir1,1).*b(ir1-1);
end %for

%% Illustrate back substitution on B using provided Matlab function
x(nref)=b(nref)/Atom(nref,2);
for ir1=nref-1:-1:1
    x(ir1)=(b(ir1)-Atom(ir1,3)).*x(ir1+1)./Atom(ir1,2);
end %for

disp('Elimination/back sub solution:  ');
disp(x);

% disp('Matlab,GNU/Octave built-in solution:  ');
% disp(A\b);
end %function 