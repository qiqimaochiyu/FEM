function data=maxwell_example3
% Yu wei 2017-05
% smooth function is known for maxwell equation
% curl(curl\textbf{u}) - \lambda\textbf{u} = \textbf{f} in \Omega
% div \textbf{u} = g  in \Omega
% t \cdot \textbf{u} = \fai on \Gamma

global w


data.f1=@f1;data.f2=@f2;
data.g=@g;data.fai=@fai;
data.exact_u1=@exact_u1;data.exact_u1x=@exact_u1x;data.exact_u1y=@exact_u1y;
data.exact_u2=@exact_u2;data.exact_u2x=@exact_u2x;data.exact_u2y=@exact_u2y;


 function r=exact_u1(x,y)
     r=-sin(2*pi*y).*(cos(2*pi*x)-1);
 end

    function r=exact_u1x(x,y)
    r=2*pi*sin(2*pi*x).*sin(2*pi*y);
    end
    function r=exact_u1y(x,y)
    r=-2*pi*cos(2*pi*y).*(cos(2*pi*x) - 1);
    end

%% U2        
    function r=exact_u2(x,y)
        r=sin(2*pi*x).*(cos(2*pi*y)-1);   
    end
    function r=exact_u2x(x,y)
        r=2*pi*cos(2*pi*x)*(cos(2*pi*y) - 1);
    end
    function r=exact_u2y(x,y)
        r=-2*pi*sin(2*pi*x)*sin(2*pi*y);
    end
%% f g boundary
    function r=g(~,~)
        r=0;
    end

    function r=fai(x,y)
        if y==0&&0<x&&x<=1
            r=exact_u1(x,y);
        elseif x==1&&0<y&&y<=1
            r=exact_u2(x,y);
        elseif y==1 && 0<=x && x<1 
            r=-exact_u1(x,y);
        elseif x==0 && 0<=y && y<1
            r=-exact_u2(x,y);
        end
    end

    function r=f1(x,y)
        r=-w*exact_u1(x,y);
    end
    function r=f2(x,y)
        r=-w*exact_u2(x,y);
    end
end


