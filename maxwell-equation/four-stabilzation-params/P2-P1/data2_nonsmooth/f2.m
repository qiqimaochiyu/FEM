function r0=f2(x,y)
global w
r=sqrt(x.^2+y.^2);
t=atan2(y,x);
t=(y>0).*t+(y<0).*(2*pi+t)+(y==0).*t;
r0=-w*2/3*r.^(-1/3).*cos(t/3);