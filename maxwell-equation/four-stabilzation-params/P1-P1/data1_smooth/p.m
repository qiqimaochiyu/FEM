function r=p(x,y)

% if -1<x && x<=0 && y==-1 || y==1 && x>=-1 && x<1
%     r=2*pi*cos(pi*x);
% elseif x==0 && -1<y && y<=0
%     r=-2*pi*cos(pi*y);
% elseif y==0 && 0<x && x<=1 
%     r=-2*pi*cos(pi*x);
% elseif x==1 && 0<y && y<=1 || x==-1&&-1<=y && y<1
%     r=2*pi*cos(pi*y);

 r=-2*pi*cos(pi*x).*cos(pi*y);   
end

