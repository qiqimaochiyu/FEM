function r=X(x,y)
global domain_type

r=sqrt(x.^2+y.^2);
t=atan2(y,x);
t=(y>0).*t+(y<0).*(2*pi+t)+(y==0).*t;
if domain_type==1
    
    if -1<x && x<=1 && y==-1 
        r=-2/3*r.^(-1/3).*sin(t/3);
    elseif x==1 && -1<y && y<=1
        r=2/3*r.^(-1/3).*cos(t/3);

    elseif y==1 && x>=-1 && x<1
        r=2/3*r.^(-1/3).*sin(t/3);
    elseif x==-1&&-1<=y && y<1
        r=-2/3*r.^(-1/3).*cos(t/3);
    end
    
elseif domain_type==2
    if -1<x&&x<=0&&y==-1 || y==0&&0<x&&x<=1
        r=-2/3*r.^(-1/3).*cos(t/3);
    elseif x==0&&-1<y&&y<=0 || x==1&&0<y&&y<=1
        r=-2/3*r.^(-1/3).*sin(t/3);

    elseif y==1 && x>=-1 && x<1
        r=2/3*r.^(-1/3).*cos(t/3);
    elseif x==-1&&-1<=y && y<1
        r=2/3*r.^(-1/3).*sin(t/3);
    end
end
end
