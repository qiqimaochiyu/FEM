function r0=u1x(x,y)

r=sqrt(x.^2+y.^2);
t=atan2(y,x);
t=(y>0).*t+(y<0).*(2*pi+t)+(y==0).*t;
r0=(2*cos(t/3).*sin(t))./(9*r^(4/3)) + (2*sin(t/3).*cos(t))./(9*r^(4/3));
