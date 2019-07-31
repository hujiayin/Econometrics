function [x, fval_neg, exitflag, output, grad, h] = maximize(fm, theta0)

%   f: the function will be maximazed

fun_neg = @(theta) -fm(theta);
%options = optimset('LargeScale','on','GradObj','on','Hessian','on');
[x, fval, exitflag, output, grad, h] = fminunc(fun_neg, theta0);
%[x, fval, exitflag, output] = fminunc(fun_neg, theta0, options);
fval_neg = -fval;

end

