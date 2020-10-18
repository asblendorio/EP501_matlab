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
    Aprime(i,1)=Atom(i,i-1);
    Aprime(i,2)=Atom(i,i);
    Aprime(i,3)=Atom(i,i+1);
end %for    
Aprime(1,2)=Atom(1,1);
Aprime(1,3)=Atom(1,2);
Aprime(50,1)=Atom(50,49);
Aprime(50,2)=Atom(50,50);
% disp('This is my Aprime matrix');
% disp(Aprime);

for ir1=2:nref
    em=Aprime(ir1,1)./Aprime(ir1-1,2);
    Aprime(ir1,1)=em;
    Aprime(ir1,2)=Aprime(ir1,2)-em.*Aprime(ir1-1,3);
    b(ir1)=b(ir1)-Aprime(ir1,1).*b(ir1-1);
end %for

%% Illustrate back substitution on B using provided Matlab function
x(nref)=b(nref)/Aprime(nref,2);
for ir1=nref-1:-1:1
    x(ir1)=(b(ir1)-Aprime(ir1,3)).*x(ir1+1)./Aprime(ir1,2);
end %for

for i=1:length(x)
    fprintf('\nx%d = %e\n',i,x(i));
end %for

finalarray = zeros(50,1);

disp('Elimination/back sub solution:  ');
disp(finalarray);

% disp('Matlab,GNU/Octave built-in solution:  ');
% disp(A\b);
end %function 