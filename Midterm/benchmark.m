%%Problem #1b
%%Produce a new version of the benchmark.(py,m) script that runs timed tests for all three of:
%%Gaussian elimination, Jacobi iteration, and your new tridiag solver. Use the same tridiagonal
%%problem for testing each of the solvers (as in the existing benchmark.(py,m) script). Store the
%%solve times for each method and each system size in separate arrays and plot times vs. number
%%of unknowns for each method on a single graph.

% Evaluate performance and scaling of Gaussian elimination, Jacobi iteration,
% and Tri-Diagonal Solver by solving systems of different size and timing the solves

nvals=50:50:500;
testtimes=zeros(size(nvals));
lrep=10;     %how many times to repeat each test

disp('Start of tests of Gaussian-elimination scaling');
for in=1:numel(nvals)
    nlarge=nvals(in);
    Blarge=diag(-1*ones(nlarge-1,1),-1)+diag(-1*ones(nlarge-1,1),1)+diag(4*ones(nlarge,1),0);    %this must be diagonally dominant or else the method won't converge
    blarge=ones(nlarge,1);
    
    for irep=1:lrep     %benchmark will repeat the same solution several times to eliminate random variations from CPU load, etc.
        tstart=cputime;
        [Blargemod,ordlarge]=Gauss_elim(Blarge,blarge);
        xlarge=backsub(Blargemod(ordlarge,:));
        tend=cputime;
        testtimes(in)=testtimes(in)+(tend-tstart)/lrep;
    end %for
    disp([' GE solution for system of size ',num2str(nlarge),' takes average time ',num2str(testtimes(in)),' s']);
end %for

figure(1);
plot(nvals,testtimes,'o','LineWidth',2,'MarkerSize',20,'MarkerFaceColor','blue')
xlabel('system size');
ylabel('time to solve (s)');
title('Empirically Determined Performance');

disp('Start of tests for Jacobi iteration');
tol=1e-9;
testtimes=zeros(size(nvals));
for in=1:numel(nvals)
    nlarge=nvals(in);
    Blarge=diag(-1*ones(nlarge-1,1),-1)+diag(-1*ones(nlarge-1,1),1)+diag(4*ones(nlarge,1),0);    %this must be diagonally dominant or else the method won't converge
    blarge=ones(nlarge,1);

    for irep=1:lrep     %benchmark will repeat the same solution several times to eliminate random variations from CPU load, etc.
        tstart=cputime;
        x0=randn(nlarge,1);
        [xit,iterations]=Jacobi(x0,Blarge,blarge,tol,false);
        tend=cputime;
        testtimes(in)=testtimes(in)+(tend-tstart)/lrep;
    end %for
    disp([' JI solution for system of size ',num2str(nlarge),' takes average time ',num2str(testtimes(in)),' s']);
end %for

figure(1);
hold on
plot(nvals,testtimes,'^','LineWidth',2,'MarkerSize',20,'MarkerFaceColor','blue')
xlabel('system size');
ylabel('time to solve (s)');
legend('Gauss elim.','Jacobi it.')
title('Empirically Determined Performance');

disp('Start of tests for Tridiagonal');
for in=1:numel(nvals)
    nlarge=nvals(in);
    Blarge=diag(-1*ones(nlarge-1,1),-1)+diag(-1*ones(nlarge-1,1),1)+diag(4*ones(nlarge,1),0);    %this must be diagonally dominant or else the method won't converge
    blarge=ones(nlarge,1);
    for irep=1:lrep     %benchmark will repeat the same solution several times to eliminate random variations from CPU load, etc.
        tstart=cputime;
        [Awork]=tridiag(Blarge,blarge);
        xlarge=backsub(Awork(1,:));
        tend=cputime;
        testtimes(in)=testtimes(in)+(tend-tstart)/lrep;
    end %for
    disp([' Tridiagonal solution for system of size ',num2str(nlarge),' takes average time ',num2str(testtimes(in)),' s']);
end %for

figure(1);
hold on
plot(nvals,testtimes,'*','LineWidth',2,'MarkerSize',20,'MarkerFaceColor','blue')
xlabel('system size');
ylabel('time to solve (s)');
legend('Gauss elim.','Jacobi it.','TriDiag')
title('Empirically Determined Performance');



