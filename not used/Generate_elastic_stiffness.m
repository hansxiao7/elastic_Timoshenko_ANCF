% Script to generate code for the evaluation of the elemental internal
% force matrix

syms L real;
syms xi yi;

syms E v;

lamda = E*v/((1+v)*(1-2*v));
miu = E/(2*(1+v));

syms xi yi;
S1 = 1 - 3*xi + 2*xi^2;
S2 = L * yi * (1 - 3*xi + 2*xi^2);
S3 = 4 * xi - 4 * xi^2;
S4 = L * yi * (4 * xi - 4 * xi^2);
S5 = -xi + 2*xi^2;
S6 = L * yi * (-xi + 2*xi^2);

S = [S1*eye(2) S2*eye(2) S3*eye(2) S4*eye(2) S5*eye(2) S6*eye(2)];
S1x = diff(S(1,:), xi) * L;
S2x = diff(S(2,:), xi) * L;
S1y = diff(S(1,:), yi) * L;
S2y = diff(S(2,:), yi) * L;

e = sym('e',[12 1]);

Sa = S1x' * S1x + S2x' *S2x;
Sb = S1y'*S1y + S2y'*S2y;
Sc = S1x'*S1y + S2x'*S2y;

syms b h;
ymax = h/(2*L);

k1 = 2*(lamda+miu)*0.5*(Sa*(e'*Sa*e-1)+Sb*(e'*Sb*e-1));
k3 = 2*miu*0.25*((Sc + Sc')*(e'*Sc*e));
KK = int(b*(k1+k3), xi, 0, 1);

K = L^2* int(KK, yi, -ymax, ymax);

matlabFunction(K, 'file', 'elastic_stiffness.m');
