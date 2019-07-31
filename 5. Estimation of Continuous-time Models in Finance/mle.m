function [thetam, fval, sigma_hat_m, p_value] = mle(fm, Xm, Ym, theta0m)

LLm = @(theta) sum(log(fm(Xm, Ym, theta)));
[thetam, fval, ~, ~, ~, h] = maximize(LLm, theta0m);

h_inv = inv(h);
sigma_hat_m = sqrt(diag(h_inv));

%ci = CI(thetam, sigma_hat_m);

p_value = PValue(thetam, sigma_hat_m);

%celldisp(ci)

%celldisp(p_value)

end

