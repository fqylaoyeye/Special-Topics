options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
fun = @example_fun;
x0 = [-1,0];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [-2,-2];
ub = [2,2];
nonlcon = @unitdisk;
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);



function [f,g] = example_fun(x)

f = 100*(x(2) - x(1)^2)^2 + (1-x(1))^2;
g = [-400*(x(2)-x(1)^2)*x(1)-2*(1-x(1));
    200*(x(2)-x(1)^2)];

end

function [c,ceq] = unitdisk(x)

c = x(1)^2 + x(2)^2 - 1;
ceq = [];

end

