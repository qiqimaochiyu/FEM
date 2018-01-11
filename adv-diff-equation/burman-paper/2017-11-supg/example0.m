function data=example0

%例子2中的数据
data.g=@g;data.a=@a;data.f=@f;
function r=g(x,y)
    if (x==0 && y<1)|| (y==0 && 0<x && x<1/2)
        r=1;
    else
        r=0;
    end
end

    function r=f(x,y)
        r=0;
    end
    function r=a(x,y)
        r=1;
    end
end