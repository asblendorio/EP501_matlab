function [invA] = LUinv(L,U)
[nref1,nref2] = size(L);
b=zeros(nref2,nref2);
invA=zeros(nref2,nref2);

k = eye(nref2);

for ir1=1:nref1
    b(1,ir1)=k(1,ir1)/L(1,1);
    for ir2=2:nref1
        j=0;
        for ir3=ir2-1:-1:1
            j=j+L(ir2,ir3)*b(ir3,ir1);
        end %for
     b(ir2,ir1)=(k(ir2,ir1)-j)/L(ir2,ir2);
    end %for
end %for

for ir1=1:nref1
    invA(nref1,ir1)=k(nref1,ir1)/U(nref1,nref1);
    for ir2=nref1-1:-1:1
        j=0;
        for ir3=ir2+1:nref1
            j=j+U(ir2,ir3)*invA(ir3,ir1);
        end %for
        invA(ir2,ir1)=b(ir2,ir1)-j/U(ir2,ir2);
    end %for
end %for
   
end %function 