function r= function_f(x,y)
if (0<=x && x<=1/2 && 0<=y && y<=1/4) || (0<=x && x<=1/4 && 0<=y && y<=1/2)
    r=32;
else
    r=0;
end
% r=logical(0<=x && x<=1/2 && 0<=y && y<=1/4) || (0<=x && x<=1/4 && 0<=y && y<=1/2);

% 
% if (1<=i && i<=9 && 1<=j && j<=17)||(1<=i && i<=17 && 1<=j && j<=9)
%     r=32;
% else
%     r=0;
% end