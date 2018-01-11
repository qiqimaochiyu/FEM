function data=example_h

%例子2中的数据
data.g=@g;data.a=@a;data.f=@f;
    function r=g(x)
        if x==3
            r=0;     
        
        elseif x>4
            r=1;    

        end
    end
    
    function r=f(~,~)
        r=0;
    end
    function r=a(~,~)
        r=1;
    end
end