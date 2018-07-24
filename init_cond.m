function [e0, ed0] = init_cond(params)
%INIT_COND Specify initial conditions for the beam nodal coordinates

% we always initialize the beam horizontal with its left end at the origin
% and at rest

e0 = zeros(params.n, 1);
ed0 = zeros(params.n, 1);

x = params.x;

for i = 1:params.ne
    e0(8*i-3, 1) = sum(x(1:i-1)) + x(i) / 2;  
    e0(8*i+1, 1) = sum(x(1:i-1)) + x(i);
end

e0(4) = 1;

for j = 1:params.ne
    e0(8*j) = 1;
    e0(8*j + 4) = 1;
end

end