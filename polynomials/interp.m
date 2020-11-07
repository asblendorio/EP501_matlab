%% Illustration of bilinear interpolation, single interval of interest

[xi,xi1,yi,yi1] = bilinear(1,1,yg,ygi,xg,xgi)
x=[xi,xi1];
y=[yi,yi1];
f=f2D;
x1=xg;
y1=yg;
% Manually written
n = length(f(:,1));
nref = length(f(:,1));

for ir1 = 1:1:n-1
    for ir2=1:1:nref-1
        [X,Y] = meshgrid(x,y);
        fvec=f(ir1:ir1+1,ir2:ir2+1);
        xvec=X(:);
        yvec=Y(:);
        M=[ones(4,1),xvec(:),yvec(:),xvec(:).*yvec(:)];
        [Mmod,order]=Gauss_elim(M,fvec(:));
        avec=backsub(Mmod(order,:));
        finterpmanual(ir1,ir2)=avec(1)+avec(2)*x1(ir1)+avec(3)*y1(ir2)+avec(4)*x1(ir1)*y1(ir2);
    end %for
end %for
k = finterpmanual;

% Matlab version
finterp=interp2(X,Y,f,x1,y1);
disp('Matlab,GNU/Octave built-in solution:');
disp(finterp);

vq = interp2(xg,yg,f2D,xgi,ygi);
plot(vq);
hold off;

%Visual interpolation point:  x1,y1
figure(1);
imagesc(x,y,f);
axis xy;
xlabel('x');
ylabel('y');
c=colorbar;
ylabel(c,'f(x,y)')
hold on;
plot(xvec,yvec,'w^','MarkerSize',15,'MarkerFaceColor','white');
plot(x1,y1,'wo','MarkerSize',20,'MarkerFaceColor','white');
hold off;

%% Illustrate cubic spline approximations using Matlab functions
x=linspace(-5,5,15);
y=sin(x);
figure(2);
plot(x,y,'o','MarkerSize',20);

splinedef=spline(x,y);
x2=linspace(min(x),max(x),256);
y2=ppval(splinedef,x2);
hold on;
plot(x2,y2,'.');

y2true=sin(x2);
plot(x2,y2true,'-');

