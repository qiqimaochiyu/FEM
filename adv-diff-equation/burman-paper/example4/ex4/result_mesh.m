global p
close all
x1=p(1,:);
x2=p(2,:);
xlin = linspace(-3,9,2401);
ylin = linspace(-3,3,1201);
[X,Y]=meshgrid(xlin,ylin);
s = X.^2+Y.^2;
flag = find( s<1 );

Z = griddata(x1,x2,uh',X,Y,'cubic');
Z(flag)=nan;

colormap('jet');
mesh(X,Y,Z);
colorbar;
figure
y11=Z(801,:);
plot(xlin,y11,'b');

figure
yl=linspace(0,3,601);
x11=Z(601:1201,1301)';
plot(yl,x11,'b');