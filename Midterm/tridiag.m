%%Problem #1a
%%Write a MATLAB or Python function tridiag() that solves a tridiagonal system of equations
%%using the Thomas algorithm. Verify your solution by applying it to the iterative test problem for
%%HW2 in the EP501 Assignments repository

function [Atom] = tridiag(A,b)
%forward elimination
nref=length(b); %system size for reference problem
Atom=cat(2,A,b); 
for ir1=2:nref
    em=Atom(ir1,1)./Atom(ir1-1,2);
    Atom(ir1,1)=em;
    Atom(ir1,2)=Atom(ir1,2)-em.*Atom(ir1-1,3);
    b(ir1)=b(ir1)-Atom(ir1,1).*b(ir1-1);
end %for
disp(Atom);

Anew = zeros(nref,nref);

Anew(1:1+N:N*N) = a;
Anew(N+1:1+N:N*N) = b;
Anew(2:1+N:N*N-N) = c;

for i=2:nref-1
    Anew(i,nref)=Atom(i,i-1);
    Anew(i,nref)=Atom(i,i);
    Anew(i,nref)=Atom(i,i+1);
end %for    
Anew(1,2)=Atom(1,1);
Anew(1,3)=Atom(1,2);
disp('New matrix for backsub');
disp(Anew);


%% Illustrate back substitution on B using provided Matlab function
% x(nref)=b(nref)/Atom(nref,2);
% for ir1=nref-1:-1:1
%     x(ir1)=(b(ir1)-Atom(ir1,3)).*x(ir1+1)./Atom(ir1,2);
% end %for
 
% disp('Dr. Z Elimination/back sub solution:  ');
% xsoln=backsub(Atom);
% disp(xsoln);

% disp('Matlab,GNU/Octave built-in solution:  ');
% disp(A\b);
end %function 