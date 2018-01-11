function r0=f1(x,y)
global w
r=sqrt(x.^2+y.^2);
t=atan2(y,x);
t=(y>0).*t+(y<0).*(2*pi+t)+(y==0).*t;
r0=-w^2*(-1/2*r.^(-1/2).*sin(t/2).*(x+1).*(y+1));
end

