function [x, fval, exitflag, output, grad, h] = maximize(fun, theta0)
%
%   fun: the function will be maximazed

fun_neg = @(theta) -fun(theta);
%options = optimset('LargeScale','on','GradObj','on','Hessian','on');
[x, fval, exitflag, output, grad, h] = fminunc(fun_neg, theta0);
%[x, fval, exitflag, output] = fminunc(fun_neg, theta0, options);

end

