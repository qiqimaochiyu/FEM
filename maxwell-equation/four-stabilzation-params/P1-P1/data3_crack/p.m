function r0=p(x,y)

r=sqrt(x.^2+y.^2);
t=atan2(y,x);
t=(y>0).*t+(y<0).*(2*pi+t)+(y==0).*t;
r0=0;