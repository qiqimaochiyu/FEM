function data=example_time

global T
%time 例子中的数据
data.g=@g;data.a=@a;data.f=@f;


    function r=a(~,~)
        r=1;
    end
    function r=f(x,y)
        if (0<=x && x<=1/2 && 0<=y && y<=1/4) || (0<=x && x<=1/4 && 0<=y && y<=1/2)
            r=1/T;
        else
            r=0;
        end
    end
    function r=g(x,y)

    if (0<=x && x<1/2 && y==0)  || (x==0 && 0<=y && y<1/2)
        r=1;
    else
        r=0;
    end
    end

end