global h_1
close all
x=linspace(0,1,1/h_1(1)+1);
y=x;
[X,Y]=meshgrid(x,y);
Z=reshape(uh,1/h_1(1)+1,1/h_1(1)+1);
figure;colormap(jet);mesh(X,Y,Z);title('new-II');view([45,15])

figure;colormap(jet);contour(X,Y,Z,20);title('new-II');colorbar;

figure;colormap(jet);B=Z(:,1/(2*h_1(1))+1)';plot(y,B,'b');title('cut line at x=0.5');

figure;colormap(jet);C=Z(1/(2*h_1(1))+1,:);plot(x,C,'b');title('cut line at y=0.5');