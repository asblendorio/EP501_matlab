%% Problem 2 (b)
%%Suppose you have a polynomial of known or given order and need to find all of its roots.
%%Write a block of code or a script that uses the exact Newton method 

i=1;
syms x;
f=@(x)(x.^5)-(15.*x.^4)+(85.*x.^3)-(225.*x.^2)+(274.*x)-120;
fprime=@(x)(5.*x.^4)-(60.*x.^3)+(255.*x.^2)-(450.*x)+274;
maxit=10;
tol=1e-6;
verbose=false;
for j=0.5:1:5.0
    x0=1;
    [root,it,success]=newton_exact(f,fprime,x0,maxit,tol,verbose);
    finalarray(1,i)=root;
    finalarray(2,i)=it;
    i=i+1;
end %for

disp(finalarray);
fplot(f,[0 6]);



