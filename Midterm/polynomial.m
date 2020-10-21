%%Problem 3b
% Write a polynomial division algorithm capable of dividing a given polynomial 
% Pn(x) (of order n and defined by a set of coefficients) by a given divisor (x − N).
% I.e. find Qn−1(x) such that: Pn(x) = (x − N)Qn−1(x) + R 
% Qn−1 is the polynomial left once you divide (x − N) out of Pn. Test your code by using it to
% divide out a factor of (x − 5) from the polynomial:
% x^5 −15x^4 +85x^3 −225x^2 +274x−120 = 0
function [b] = polynomial(a,x) 
% p = x^5-15*x^4+85*x^3-225*x^2+274*x-120;
% f = x - 5;
% a = [1 -15 85 -225 274 -120];
% x = [1 -5];
% q = polynomialReduce(coef,f);
% p(x)=a(1) x^n +a(2) x^(n-1) + ...
nref = length(a);
b=zeros(1,nref-1);
b(1)=a(1);
for i=1:nref-1
    b(i+1)= x.*b(i) + a(i+1);
end %for  

%f=(z.^(nref-1)).*f;

% b=zeros(1,nref-1); 
% for i=1:nref-2
%     b(i+1)=z.*b(i)+a(i+1);
% end %for 
disp('The polynomial is reduced to:');
disp(b);
end %function  
 
