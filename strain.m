function epislon = strain(e, x, y, L)
% This function is to calculate the strains in the shear-deformable
% ANCF beam given positions (x,y). The returned strain vector has three
% components strain = [strain_xx, strain_yy, strain_xy];

Sx = shape_fun(x, y, L, 1);
Sy = shape_fun(x, y, L, 2);

Sx = [Sx(1)*eye(2), Sx(2)*eye(2),Sx(3)*eye(2),Sx(4)*eye(2),Sx(5)*eye(2), Sx(6)*eye(2)];
Sy = [Sy(1)*eye(2), Sy(2)*eye(2),Sy(3)*eye(2),Sy(4)*eye(2),Sy(5)*eye(2), Sy(6)*eye(2)];
S1x = Sx(1, :);
S2x = Sx(2,:);
S1y = Sy(1,:);
S2y = Sy(2,:);

Sa = S1x' * S1x + S2x' * S2x;
Sb = S1y' * S1y + S2y' * S2y;
Sc = S1x' * S1y + S2x' * S2y;

strain_1 = 0.5 * (e'*Sa*e - 1);
strain_2 = 0.5 * (e'*Sb*e - 1);
strain_3 = 0.5 * e' * Sc *e;

epislon = [strain_1; strain_2; strain_3];

end