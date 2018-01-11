function r0=f2(x,y)
global
r=sqrt(x.^2+y.^2);
t=atan2(y,x);
t=(y>0).*t+(y<0).*(2*pi+t)+(y==0).*t;
r0=-w^2*1/2*r.^(-1/2).*cos(t/2).*(x+1).*(y+1);
end
