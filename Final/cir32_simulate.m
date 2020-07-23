%% Function to Simulate CIR 3/2 Data
function data = cir32_simulate(a, b, c, h, size, initial_value)
    X = zeros(size,1);
    X(1) = initial_value;
    for i = 2:size
        X(i) =  b * (a- X(i-1)) *h + c * (X(i-1)^(3/2)) * sqrt(h) * normrnd(0,1) + X(i-1);
    end
    data = X;
end

