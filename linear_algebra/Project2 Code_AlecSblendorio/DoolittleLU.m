%%Problem #1 
%%Create a new version of your simple forward elimination function from the 
%%first assignment so that it performs Doolittle LU factorization.
 
function [Awork,L,bprime,x] = DoolittleLU(A,b)
nref=length(b);                %system size for reference problem
%Compute U
%note that the elimination procedure coded below modifies the matrix B
Awork=cat(2,A);
L=eye(nref,nref);
                                                       %This is our working version of the matrix used to perform elimination (i.e. it will be modified)
for ir1=2:nref                                           %loop over rows from 2 to n performing elimination, this index marks what row we are starting the elimination from (i.e. using) for this particular column
    for ir2=ir1:nref                                        %this index marks the present position where elimination is being performed - i.e. where we are applying the elementary row operations
    fact=Awork(ir2,ir1-1);
    Awork(ir2,:)=Awork(ir2,:)-fact/Awork(ir1-1,ir1-1).*Awork(ir1-1,:);
    %multiplier of the variable we are attempting to eliminate, its ir-1 column of this row    
    %store elimination multipliers below                                    %multiplier of the variable we are attempting to eliminate, its ir-1 column of this row
    fact1 = fact/Awork(ir1-1,ir1-1);                        % numerator is the previous row previous column. 
    L(ir2,ir1-1)=fact1;
    end %for 
end %for

U=Awork;

disp('Upper Triangular Forward Elim([Aref,bref]) = ');
disp(U);
disp('Lower Triangular Forward Elim([Aref,bref]) = ');
disp(L);

%LU Factorization
LU=L*U;
disp('The LU Factorization of the Matrix A is:'); disp(LU);
%% Forward Substitution Solution 
bprime=fwdsub1(L,b);
disp('Forward Substitution Results/B Prime Vector =');
disp(bprime);

%% Back substitution solution
xsoln=backsub(U,bprime);
disp('Elimination/back sub solution:  ');
disp(xsoln);


end %function 



