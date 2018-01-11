function r=q(x,y)
% Yu wei 2017-05
% L-shape 
%  1  _ _ _ _ _ _ _ _1
%    |              |
%    |              |
%  0 |       _ _ _ _| 0     
%    |       |0
%    |       |
%  -1|_ _ _ _|-1
global domain_type
if domain_type==1
    if -1<x && x<=1 && y==-1 
        r=cos(pi*y)*sin(pi*x);
    elseif x==1 && -1<y && y<=1
        r=cos(pi*x)*sin(pi*y);
    elseif y==1 && x>=-1 && x<1
        r=-cos(pi*y)*sin(pi*x);
    elseif x==-1&&-1<=y && y<1
        r=-cos(pi*x)*sin(pi*y);
    end
elseif domain_type==2
    
    if -1<x&&x<=0&&y==-1 || y==0&&0<x&&x<=1
        r=cos(pi*y)*sin(pi*x);
    elseif x==0&&-1<y&&y<=0 || x==1&&0<y&&y<=1
        r=cos(pi*x)*sin(pi*y);
    elseif y==1 && x>=-1 && x<1
        r=-cos(pi*y)*sin(pi*x);
    elseif x==-1&&-1<=y && y<1
        r=-cos(pi*x)*sin(pi*y);
    end
end

end
