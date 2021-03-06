function p = ancf_params()
%ANCF_PARAMS Problem specification, define the sections

L = 2;%total beam length;
p.ne = 20;% number of beam elements
for i = 1:p.ne
    p.x(i) = L/p.ne;   % noodle length
end

p.n = 8 * p.ne + 4;  % number of nodal coordinates
p.b = 0.1;
p.h = 0.5;
p.rho = 2.17e-7;      % density k/in^3
p.E = 2.07e11; %elastic modulus
p.v = 0.3; %Poisson's ratio
% p.I = 34*26^3/12;
p.F = 500e6* p.h^3;


% specify constraints at the noodle ends
%the element in vector represent which dof is constrained, from 1 to 4
p.ml = [1,2,3,4]; %left end
p.mr = [];%right end
end