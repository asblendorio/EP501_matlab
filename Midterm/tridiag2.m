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
    Aprime(i,i)=Atom(i,i);
    Aprime(i,i+1)=Atom(i,i+1);
end %for    
Aprime(1,1)=Atom(1,1);
Aprime(1,3)=Atom(1,2);
Aprime(50,49)=Atom(50,49);
Aprime(50,50)=Atom(50,50);

for ir1=2:nref
    em=Aprime(ir1,1)./Aprime(ir1-1,2);
    Aprime(ir1,1)=em;
    Aprime(ir1,2)=Aprime(ir1,2)-em.*Aprime(ir1-1,3);
    b(ir1)=b(ir1)-Aprime(ir1,1).*b(ir1-1);
end %for

disp('This is my Aprime matrix');
disp(Aprime);

%% Illustrate back substitution on B using provided Matlab function

xsoln=backsub(Aprime);
disp('Elimination/back sub solution:  ');
disp(xsoln);

x(nref)=b(nref)/Aprime(nref,2);
for ir1=nref-1:-1:1
    x(ir1)=(b(ir1)-Aprime(ir1,3)).*x(ir1+1)./Aprime(ir1,2);
end %for
% 
% % for i=1:length(x)
% %     fprintf('\nx%d = %e\n',i,x(i));
% % end %for
% x = x.';
% disp('Elimination/back sub solution:  ');
% fprintf('%e\n',x);

% disp('Matlab,GNU/Octave built-in solution:  ');
% disp(a\b);
end %function 
