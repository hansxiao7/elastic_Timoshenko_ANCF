clc
clear
%This function is to calculate the static response with the Newton-Raphson
%method
params = ancf_params();

F = zeros(params.n, 1);
F(params.n-2, 1) = -params.F;

[e, ~] = init_cond(params);

for j = 1:10000
    j
    Qint = total_internal_force(e, params);
    Q_gradient= Q_gradient_total(e, params);
    K = Q_gradient(5:params.n, 5:params.n);
    Kinv = inv(K);
    
    delta_e = Kinv*(F(5:params.n) - Qint(5:params.n));
    
    e = e + [0;0;0;0;delta_e];
    
end