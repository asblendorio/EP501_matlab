%% Integration using the Iterated Trapezoidal method
lx=50;
ly=50;
lz=50;

%constants
Q=1;
a=1;
e0=8.854.*10.^(-12);
phi=zeros(100,100,100);

x = linspace(-3*a,3*a,lx);
y = linspace(-3*a,3*a,ly);
z = linspace(-3*a,3*a,lz);
[X,Y,Z]=meshgrid(x,y,z);

% dp_dx=zeros(size(phi));
% dp_dy=zeros(size(phi));
% dp_dz=zeros(size(phi));

%gradient of scalar function
dx=x(2)-x(1);
dy=y(2)-y(1);
dz=z(2)-z(1);
DE = 1; % energy initiliaztion 

for i=1:50
    for j=1:50
        for k=1:50
            phi(i,j,k) = Q/(4.*pi.*e0.*sqrt(x(i).^2 +y(j).^2 + z(k).^2)); %boundary conditions sets this value at zero outside of source field
            % Integral 
            DE = ((0.5.*phi(i+1,j,k)+phi(i,j,k).*dx)+ ...
            (0.5.*phi(i,j+1,k)+phi(i,j,k).*dy)+...
            (0.5.*phi(i,j,k+1)+phi(i,j,k).*dz) + DE);
           
        end %for  
    end %for
end %for
We = DE;
fprintf('\n Total electrostatic energy We = %e Joules\n',We);

