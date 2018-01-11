function r=exact_solution(x,y)


r=x.*y.*(1-x/2).*(1-y).*exp(x+y);