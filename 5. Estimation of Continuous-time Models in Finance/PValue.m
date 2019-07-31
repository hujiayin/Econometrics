function [p_value] = PValue(thetam, sigma_hat_m)

[m, n] = size(thetam);
for i = 1:1:n
    p_value(i) = 2 * (1-normcdf(abs(thetam(i)), 0, sigma_hat_m(i)));
end

end

