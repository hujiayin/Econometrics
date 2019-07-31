function mle(fm, Xm, Ym, theta0m, theta0)

LLm = @(theta) sum(log(fm(Xm, Ym, theta)));
[thetam, fval, exitflag, output, grad, h] = maximize(LLm, theta0m);

h_inv = inv(h);
sigma_hat_m = sqrt(diag(h_inv));

ci = CI(thetam, sigma_hat_m);

p_value = PValue(thetam, theta0, sigma_hat_m);

thetam

fval

celldisp(ci)

celldisp(p_value)

end

