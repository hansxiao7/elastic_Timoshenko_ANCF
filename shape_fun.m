function S = shape_fun(x, y, L, der)
%SHAPE_FUN calculate the shape functions at the specified location
%   This function returns the shape functions or their derivatives
%   for a shear-deformable ANCF beam element of length L
%   at the specified curvilinear coordinate x along the element axis
%   (0 <= x <= L) and the coodinate y at the cross-section (-h/2 <= y <= h/2).
%   S = SHAPE_FUN(x, y, L, 0)  returns S
%   S = SHAPE_FUN(x, y, L, 1)  returns Sx
%   S = SHAPE_FUN(x, y, L, 2)  returns Sy

xi = x / L;
yi = y / L;

switch der
    case 0
        S = [1 - 3*xi + 2*xi^2
            L * yi * (1 - 3*xi + 2*xi^2)
            4* xi - 4 * xi^2
            L * yi * (4* xi - 4 * xi^2)
            -xi + 2* xi^2
            L * yi * (-xi + 2* xi^2)];
    case 1
        S = [(-3 + 4 * xi)/L
            yi * (-3 + 4 * xi)
            (4 - 8 * xi) / L
            yi * (4 - 8 * xi)
            (-1 + 4 * xi)/ L
            yi * (-1 + 4 * xi)];
    case 2
        S = [0
             1 - 3*xi + 2*xi^2
             0
             4* xi - 4 * xi^2
             0
             -xi + 2* xi^2];
end

end