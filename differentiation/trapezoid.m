%% Integration in multiple dimensions
%using an iterated trapezoidal method (sweeps of single dimensional integrations) or multi-dimensional
%trapezoidal method program that you write.
lx=100;
ly=100;
lz=100;

%constants
Q=1;
a=1;
e0=8.854.*10.^(-12);

x = linspace(-3*a,3*a,lx);
y = linspace(-3*a,3*a,ly);
z = linspace(-3*a,3*a,lz);
[X,Y,Z]=meshgrid(x,y,z);

phi=zeros(size(X));
%terms inside integral
int=(e_0*laplace_phi).*phi;
s=zeros(size(Y,1));


%Integrating F at every y and z for all values of x 
for k=1:size(Z,1)-1
    for j=1:size(Y,3)-1 
        for i=2:size(X,2)-1
        s(j,k)=s(j,k)+(2*int(j,i,k));
        end
        s(j,k)=s(j,k)+int(j,1,k)+int(j,size(X,2),k);
    end
end

plane=zeros(size(Z,3),1);
%Integration for values of y at every z 
for k=1:size(Z,1)
    for j=2:size(Y,3)-1
    plane(k)=plane(k)+(2*surface(j,k));
    end
    plane(k)=plane(k)+surface(1,k)+surface(size(Y,3),k);             
    %Integration for last and first term
end

%Integrating for all values of z 
for k=2:size(Z,1)-1
    Integrated_term=Integrated_term+(2*plane(k));
end

Integrated_term=Integrated_term+plane(1)+plane(size(Z,3)); We=(-0.5)*Integrated_term;
fprintf('\n Total electrostatic energy We = %f Joules\n',We);
Total electrostatic energy We = 4978705645790294.000000 Joules
