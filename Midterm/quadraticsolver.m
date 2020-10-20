%%use quadratic function to solve the column vector of coefficients [a, b, c] .
%%Test your code on the quadratic:
%%2x^2 − 6x+4=0 
syms x
coef=[2;-6;4];
f=@quadratic;
x=f(coef,x);
disp(x);


