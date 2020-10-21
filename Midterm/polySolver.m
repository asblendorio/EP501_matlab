% Problem 3d and e
% 
function[r,Pnew]= polySolver(a,r0,tol,nmax)
nref=length(a)-1;
r=zeros(nref,1);
Pnew=zeros(nref,1);
% r0=0;
% tol=1e-9;
% nmax=100;
for i=1:1:nref
    ni=0;
    x=r0;
    err=10;
    while ni<nmax && err>=tol 
        [v,p]=polynomial(a,x);
        [dv,~]=polynomial(p,x);
        xn=x-v/dv;
        err=abs(xn-x);
        ni=ni+1;
        x=xn;
    end %while    
    r(i)=x;
    Pnew(i)=ni;
    [~,a]=polynomial(a,x); %deflated
      
end %for

end %function 