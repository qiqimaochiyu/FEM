function r0=q2(x,y)

r=sqrt(x.^2+y.^2);
t=atan2(y,x);
t=(y>0).*t+(y<0).*(2*pi+t)+(y==0).*t;

if x==0&&y==0
    r0=1/2*r.^(-1/2).*cos(t/2).*(x+1).*(y+1);
elseif y==0&&0<x&&x<=1
    r0=-1/2*r.^(-1/2).*cos(t/2).*(x+1).*(y+1);
end

 