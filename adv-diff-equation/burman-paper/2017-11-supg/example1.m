function data=example1

%����1�е�����
data.g=@g;data.a=@a;data.f=@f;
function r=g(x,y)
    if (y==0 && x>0)|| (x==1 && 0<y && y<(1/2*tan(55/180*pi)))
        r=0;
    else
        r=1;
    end
end

    function r=f(x,y)
        r=0;
    end
    function r=a(x,y)
        r=1;
    end
end
