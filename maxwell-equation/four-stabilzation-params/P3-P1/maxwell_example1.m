function data=maxwell_example1
% Yu wei 2017-05
% smooth function is known for maxwell equation
% curl(curl\textbf{u}) - \omega^2\textbf{u} = \textbf{f} in \Omega
% div \textbf{u} = g  in \Omega
% curl \textbf{u} = \lambda, \textbf{u} \cdot \textbf{n} = \fai on \Gamma


data.omega=0;
data.f1=@f1;data.f2=@f2;
data.g=@g;data.lambda=@lambda;data.fai=@fai;
data.exact_u1=@exact_u1;data.exact_u1x=@exact_u1x;data.exact_u1y=@exact_u1y;
data.exact_u2=@exact_u2;data.exact_u2x=@exact_u2x;data.exact_u2y=@exact_u2y;
%% U1
    function r=exact_u1(x,y)
        r=sin(pi*y).*cos(pi*x);
    end
    function r=exact_u1x(x,y)
        r=-pi*sin(pi*x).*sin(pi*y);
    end
    function r=exact_u1y(x,y)
        r=pi*cos(pi*x).*cos(pi*y);
    end

%% U2        
    function r=exact_u2(x,y)
        r=-sin(pi*x).*cos(pi*y);
    end
    function r=exact_u2x(x,y)
        r=-pi*cos(pi*x).*cos(pi*y);
    end
    function r=exact_u2y(x,y)
        r=pi*sin(pi*x).*sin(pi*y);
    end
%% f g boundary
    function r=g(~,~)
        r=0;
    end

    function r=lambda(x,y)
        r=exact_u2x(x,y)-exact_u1y(x,y);
    end

    function r=fai(x,y)
        if y==-1 && -1<x && x<=0 || y==0 && 0<x && x<=1
            r=-exact_u2(x,y);
        elseif x==0 && -1<y && y<=0 || x==1 && 0<y && y<=1
            r=exact_u1(x,y);
        elseif y==1 && -1<=x && x<1 
            r=exact_u2(x,y);
        elseif x==-1 && -1<=y && y<1
            r=-exact_u1(x,y);
        end
    end

    function r=f1(x,y)
        r=2*pi^2*cos(pi*x).*sin(pi*y)-data.omega^2*exact_u1(x,y);
    end
    function r=f2(x,y)
        r=-2*pi^2*cos(pi*y)*sin(pi*x)-data.omega^2*exact_u1(x,y);
    end
        
end

