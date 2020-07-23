%% Function to Simulate square Radial OU Data
function data = ou_simulate(beta, sigma, h, size, initial_value)
    X = zeros(size,1);
    X(1) = initial_value;
    for i = 2:size
        X(i) = 1 + 2 * beta * X(i-1) * h + 2 * sigma * sqrt(X(i-1) * h)* normrnd(0,1) + X(i-1);
    end
    data = X;
end