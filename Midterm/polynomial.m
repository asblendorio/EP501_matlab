%%Problem 3b
% Write a polynomial division algorithm capable of dividing a given polynomial 
% Pn(x) (of order n and defined by a set of coefficients) by a given divisor (x − N).
% I.e. find Qn−1(x) such that: Pn(x) = (x − N)Qn−1(x) + R 
% Qn−1 is the polynomial left once you divide (x − N) out of Pn. Test your code by using it to
% divide out a factor of (x − 5) from the polynomial:
% x^5 −15x^4 +85x^3 −225x^2 +274x−120 = 0
function [q] = polynomial(coef,f)
syms x n
% p = x^5-15*x^4+85*x^3-225*x^2+274*x-120;
% f = x - 5;
% coef = [1 -15 85 -225 274 -120];
% f = [1 -5];
q = polynomialReduce(coef,f);

disp(q)

