%%Problem #1 
%%Create a new version of your simple forward elimination function from the 
%%first assignment so that it performs Doolittle LU factorization.
%testing new branch edit

function [Awork,Alow] = DoolittleLU(A,L,b)
nref=length(b);                %system size for reference problem
% Compute U
%note that the elimination procedure coded below modifies the matrix B
Awork=cat(2,A,b);          %This is our working version of the matrix used to perform elimination (i.e. it will be modified)
for ir1=2:nref                                           %loop over rows from 2 to n performing elimination, this index marks what row we are starting the elimination from (i.e. using) for this particular column
    for ir2=ir1:nref                                     %this index marks the present position where elimination is being performed - i.e. where we are applying the elementary row operations
        fact=Awork(ir2,ir1-1);                                    %multiplier of the variable we are attempting to eliminate, its ir-1 column of this row
        Awork(ir2,:)=Awork(ir2,:)-fact/Awork(ir1-1,ir1-1).*Awork(ir1-1,:);    %subtract off previous row modified by a factor that eliminates the ir-1 column term in this row (so it has only super-diagonal elements), this is a little bit wasteful as it uses entire row...
    end %for
end %for

disp('elim([Aref,bref]) = ');
disp(Awork);

x = zeros(nref,1);
Alow=cat(2,L,b);
x(1)=Alow(1,nref+1)/Alow(1,1);

%Compute L
for i=2:nref
    x(i,1)=(b(i)-L(i,1:i-1)*x(1:i-1,1))./L(i,i);    
end %for 

 Alow(nref) = b(nref)/L(nref,nref);
 disp('Lower Triangular Forward Elim([Aref,bref]) = ');
 disp(Alow);

 %LU Factorization
 LU=Awork.*Alow;
 disp('The LU Factorization of the Matrix A is:');
 disp(LU);
%% Back substitution solution
xsoln=backsub(LU);
disp('Elimination/back sub solution:  ');
disp(xsoln);

%% Solve the test linear system of equations using LU Factorization and Back-Sub
% d = Alow\b;
% x = Awork\d;
% disp(x);


end %function 



