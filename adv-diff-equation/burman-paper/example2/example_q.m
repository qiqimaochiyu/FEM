function data=example_q

%例子2中的数据
data.g=@g;data.a=@a;data.f=@f;
function r=g(x,y)
    if (x==0 && y<1)|| (y==0 && 0<x && x<1/2)
        r=1;
    else
        r=0;
    end
end

    function r=f(~,~)
        r=0;
    end
    function r=a(~,~)
        r=1;
    end
end