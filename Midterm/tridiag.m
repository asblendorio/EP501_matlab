%%Problem #1a
%%Write a MATLAB or Python function tridiag() that solves a tridiagonal system of equations
%%using the Thomas algorithm. Verify your solution by applying it to the iterative test problem for
%%HW2 in the EP501 Assignments repository

function [Atom,ord] = tridiag(Ait,bit)
nref=length(bit);                %system size for reference problem
Atom=cat(2,Ait,bit); 
for ir1=2:nref
    for ir2=ir1:nref
        fact=Atom(ir2,ir1-1);                                       
        Atom(ir2,:)=Atom(ir2,:)-fact/Atom(ir1-1,ir1-1).*Atom(ir1-1,:);
    end %for
end %for

disp('elim([Aref,bref]) = ');
disp(Atom);
%% Illustrate back substitution on B using provided Matlab function
xsoln=backsub(Atom);
disp('Elimination/back sub solution:  ');
disp(xsoln);
disp('Matlab,GNU/Octave built-in solution:  ');
disp(Ait\bit);
end %function 