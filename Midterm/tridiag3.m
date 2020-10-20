%Problem #1a
%%Write a MATLAB or Python function tridiag() that solves a tridiagonal system of equations
%%using the Thomas algorithm. Verify your solution by applying it to the iterative test problem for
%%HW2 in the EP501 Assignments repository

function [Jerry] = tridiag3(A,b)
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
% Aprime(50,1)=Atom(50,49);
% Aprime(50,2)=Atom(50,50);
Anew = zeros(50,50);

for ir1=2:nref-1
    em=Aprime(ir1,1)./Aprime(ir1-1,2);
    Anew(ir1,ir1)=Aprime(ir1,2)-em.*Aprime(ir1-1,3);
    Aprime(ir1,2)=Aprime(ir1,2)-em.*Aprime(ir1-1,3);
    Anew(ir1,ir1+1)=-1;
    %b(ir1)=b(ir1)-Aprime(ir1,1).*b(ir1-1);
end %for
Anew(1,1) = 4;
Anew(1,2) = -1;
disp('Aprime');
disp(Aprime);
disp('This is my new matrix');
disp(Anew);
Jerry=cat(2,Anew,b);
disp('Jerry matrix:');
disp(Jerry);
%% Illustrate back substitution on B using provided Matlab function
xback = backsub(Jerry);
disp('Back Substitution of matrix:');
disp(xback);
% x(nref)=b(nref)/Aprime(nref,2);
% for ir1=nref-1:-1:1
%     x(ir1)=(b(ir1)-Aprime(ir1,3)).*x(ir1+1)./Aprime(ir1,2);
% end %for

% for i=1:length(x)
%     fprintf('\nx%d = %e\n',i,x(i));
% end %for
% x = x.';
% disp('Elimination/back sub solution:  ');
% fprintf('%e\n',x);

% disp('Matlab,GNU/Octave built-in solution:  ');
% disp(a\b);
end %function 