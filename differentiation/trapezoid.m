%% Multidimensional function and partial derivatives:  grad and div (laplacian and curl also useful)
lx=100;
ly=100;
lz=100;

%constants
Q=1;
a=1;
e0=8.854.*10.^(-12);
phi=zeros(100,100,100);

x = linspace(-3*a,3*a,lx);
y = linspace(-3*a,3*a,ly);
z = linspace(-3*a,3*a,lz);
[X,Y,Z]=meshgrid(x,y,z);

dp_dx=zeros(size(phi));
dp_dy=zeros(size(phi));
dp_dz=zeros(size(phi));

%gradient of scalar function
dx=x(2)-x(1);
dy=y(2)-y(1);
dz=z(2)-z(1);
DE = 1; % energy initiliaztion 

for i=1:100
    for j=1:100
        for k=1:100
            if (sqrt(x(i).^2 +y(j).^2+z(k).^2) < a)
            f(i,j,k) = (-1/2.*e0.*-6.*(Q/(8.*pi.*e0.*a.^3))).*(Q/(4.*pi.*e0.*a)-(Q/(8.*pi.*e0.*a.^3))*sqrt(x(i).^2 +y(j).^2+z(k).^2).^2-a.^2);
            
            else
            f(i,j,k) = 0; %boundary conditions sets this value at zero outside of source field
            end %if
            % Integral 
            DE = (0.5.*f(i+1,j,k)+f(i,j,k).*dx)+ (0.5.*f(i,j+1,k)+f(i,j,k).*dy)+(0.5.*f(i,j,k+1)+f(i,j,k).*dz) + DE;
        end %for  
    end %for
end %for

disp(DE);

fprintf('\n Total electrostatic energy We = %f Joules\n',We);
Total electrostatic energy We = 4978705645790294.000000 Joules
