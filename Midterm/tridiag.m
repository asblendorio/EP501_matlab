%%Problem #1
%%Write a MATLAB or Python function tridiag() that solves a tridiagonal system of equations
%%using the Thomas algorithm. Verify your solution by applying it to the iterative test problem for
%%HW2 in the EP501 Assignments repository

function [Atom] = tridiag(A,b)
nref=length(b);                %system size for reference problem

Atom=cat(2,A,b);          %This is our working version of the matrix used to perform elimination (i.e. it will be modified)
for ir1=2:nref                                           %loop over rows from 2 to n performing elimination, this index marks what row we are starting the elimination from (i.e. using) for this particular column
    for ir2=ir1:nref                                     %this index marks the present position where elimination is being performed - i.e. where we are applying the elementary row operations
        fact=Atom(ir2,ir1-1);                                    %multiplier of the variable we are attempting to eliminate, its ir-1 column of this row
        Atom(ir2,:)=Atom(ir2,:)-fact/Atom(ir1-1,ir1-1).*Atom(ir1-1,:);    %subtract off previous row modified by a factor that eliminates the ir-1 column term in this row (so it has only super-diagonal elements), this is a little bit wasteful as it uses entire row...
    end %for
end %for

%disp('elim([Aref,bref]) = ');
%disp(Atom);
%% Illustrate back substitution on B using provided Matlab function
xsoln=backsub(Atom);
disp('Elimination/back sub solution:  ');
disp(xsoln);
disp('Matlab,GNU/Octave built-in solution:  ');
disp(A\b);
end %function 