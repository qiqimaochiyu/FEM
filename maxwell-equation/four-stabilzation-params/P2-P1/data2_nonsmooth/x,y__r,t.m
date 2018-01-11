
r=sqrt(x.^2+y.^2);
t=atan2(y,x);
t=(y>0).*t+(y<0).*(2*pi+t)+(y==0).*t;

[diff(u,x) diff(u,y)] = 
[diff(u,r) diff(u,t)]*
[cos(t)          sin(t);
-1/r*sin(t)  1/r*cos(t)];


f=0;g=0;p=0;