% Script to generate code for the evaluation of the elemental mass matrix
% for a shear-deformable ANCF beam element.
% Note: we assume that the density is constant over the element.

syms L h real;

syms xi yi;
S1 = 1 - 3*xi + 2*xi^2;
S2 = L * yi * (1 - 3*xi + 2*xi^2);
S3 = 4 * xi - 4 * xi^2;
S4 = L * yi * (4 * xi - 4 * xi^2);
S5 = -xi + 2*xi^2;
S6 = L * yi * (-xi + 2*xi^2);

S = [S1*eye(2) S2*eye(2) S3*eye(2) S4*eye(2) S5*eye(2) S6*eye(2)];
ymax = h/(2*L);

M = int(S'*S, xi, 0, 1);
M = simplify(int(M, yi, -ymax, ymax));

matlabFunction(M, 'file', 'mass.m');