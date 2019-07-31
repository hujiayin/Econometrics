function [ci] = CI(thetam, sigma_hat_m)

conf_level = 0.05;
[m, n] = size(thetam);
l_m = norminv(1-conf_level/2,0, sigma_hat_m);


for i = 1:1:n
    ci(i) = {[thetam(i) - l_m(i), thetam(i) + l_m(i)]};
end


end

