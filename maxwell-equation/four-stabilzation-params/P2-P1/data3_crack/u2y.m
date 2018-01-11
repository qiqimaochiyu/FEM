function r0=u2y(x,y)

r=sqrt(x.^2+y.^2);
t=atan2(y,x);
t=(y>0).*t+(y<0).*(2*pi+t)+(y==0).*t;
r0=- (cos(t/2).*sin(t).*(x + 1).*(y + 1))./(4*r.^(3/2)) - (sin(t/2).*cos(t).*(x + 1).*(y + 1))./(4*r.^(3/2));

end
