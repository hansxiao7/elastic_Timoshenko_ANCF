function Qint = internal_force(e, L, params)

C = elastic_c_matrix(params);
h = params.h;

[x5, w5] = gauss_points(5);
y5 = x5;

Qint = zeros(12, 1);

for i = 1:5
    x = (x5(i)+1)*L/2;
    for j = 1:5
        y = y5(j)*h/2;
        epislon = strain(e, x, y, L);
        Sx = shape_fun(x, y, L, 1);
        Sy = shape_fun(x, y, L, 2);
        Sx = [Sx(1)*eye(2) Sx(2)*eye(2) Sx(3)*eye(2) Sx(4)*eye(2) Sx(5)*eye(2) Sx(6)*eye(2)];
        Sy = [Sy(1)*eye(2) Sy(2)*eye(2) Sy(3)*eye(2) Sy(4)*eye(2) Sy(5)*eye(2) Sy(6)*eye(2)];
        Sx1 = Sx(1,:);
        Sx2 = Sx(2,:);
        Sy1 = Sy(1,:);
        Sy2 = Sy(2,:);
        Sa = Sx1'*Sx1 + Sx2'*Sx2;
        Sb = Sy1'*Sy1 + Sy2'*Sy2;
        Sc = Sx1'*Sy1 + Sx2'*Sy2;
        strain_gradient = [Sa*e, Sb*e, 0.5*(Sc'+Sc)*e];
        Qint = Qint + params.b * w5(i)*w5(j)*strain_gradient*C*epislon * L/2 * h/2;

    end
end

end

