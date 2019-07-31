function [p_value] = PValue(thetam, theta0, sigma_hat_m)

[m, n] = size(thetam);
for i = 1:1:n
    p_value(i) = {2 * (1-normcdf(abs(thetam(i)-theta0(i)), 0, sigma_hat_m(i)))};
end

end

