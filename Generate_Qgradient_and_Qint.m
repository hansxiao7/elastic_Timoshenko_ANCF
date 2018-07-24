%This function is to generate the gradient of the internal force vector
%Qint to use it in the static analysis

e = sym('e', [12 1], 'real');

syms E v real;
C = E/(1-v^2)*[1 v 0; v 1 0; 0 0 (1-v)/2];

syms x y real;
syms L real;

xi = x/L;
yi = y/L;

Sx =[(-3 + 4 * xi)/L
     yi * (-3 + 4 * xi)
     (4 - 8 * xi) / L
     yi * (4 - 8 * xi)
     (-1 + 4 * xi)/ L
     yi * (-1 + 4 * xi)];
 
Sy = [0
       1 - 3*xi + 2*xi^2
       0
       4* xi - 4 * xi^2
       0
       -xi + 2* xi^2];
   
Sx = [Sx(1)*eye(2) Sx(2)*eye(2) Sx(3)*eye(2) Sx(4)*eye(2) Sx(5)*eye(2) Sx(6)*eye(2)];
Sy = [Sy(1)*eye(2) Sy(2)*eye(2) Sy(3)*eye(2) Sy(4)*eye(2) Sy(5)*eye(2) Sy(6)*eye(2)];
Sx1 = Sx(1,:);
Sx2 = Sx(2,:);
Sy1 = Sy(1,:);
Sy2 = Sy(2,:);
Sa = Sx1'*Sx1 + Sx2'*Sx2;
Sb = Sy1'*Sy1 + Sy2'*Sy2;
Sc = Sx1'*Sy1 + Sx2'*Sy2;

epislon = [0.5*(e'*Sa*e-1);0.5*(e'*Sb*e-1);0.5* e'*Sc*e];

strain_gradient = [];
for i = 1:12
    strain_gradient = [strain_gradient, simplify(diff(epislon, e(i)))];
end

Qint = strain_gradient' * C * epislon;

Q_gradient = [];
for j = 1:12
    Q_gradient = [Q_gradient, simplify(diff(Qint, e(j)))];
end

matlabFunction(Qint, 'file', 'Qint.m');
matlabFunction(Q_gradient, 'file', 'Q_gradient.m');


