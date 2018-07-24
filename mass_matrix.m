function M = mass_matrix(params)
%MASS_MATRIX Calculate the generalized mass matrix for a beam.
%ml and mr are the constraints at beam end dof, for example, ml = [1,2,3,4]
%means a fixed end
% Evaluate the elemental mass matrix.
M = zeros(params.n);
ml = params.ml;
mr = params.mr;

for i = 1:max(size(params.x))
    L = params.x(i);
    Me = params.rho * params.b * L^2 *mass(L, params.h);
    
    % Assemble mass matrix for a noodle with 'ne' elements.
    M00 = Me(1:4,1:4);
    M01 = Me(1:4,5:8);
    M02 = Me(1:4,9:12);
    M10 = Me(5:8,1:4);
    M11 = Me(5:8,5:8);
    M12 = Me(5:8,9:12);
    M20 = Me(9:12,1:4);
    M21 = Me(9:12,5:8);
    M22 = Me(9:12,9:12);
    
    ii = (i-1)*8 + 1;
    M(ii:ii+3, ii:ii+3) = M(ii:ii+3, ii:ii+3) + M00;
    M(ii:ii+3, ii+4:ii+7) = M01;
    M(ii:ii+3, ii+8:ii+11) = M02;
    M(ii+4:ii+7, ii:ii+3) = M10;
    M(ii+4:ii+7, ii+4:ii+7) = M11;  
    M(ii+4:ii+7,ii+8:ii+11) = M12;
    M(ii+8:ii+11, ii:ii+3) = M20;
    M(ii+8:ii+11, ii+4:ii+7) = M21;
    M(ii+8:ii+11, ii+8:ii+11) = M(ii+8:ii+11, ii+8:ii+11) + M22;
end
% If some nodal coordinates are to be fixed, modify the mass matrix
% (zero out rows and columns, then set diagonal block(s) to identity)


if max(size(ml)) ~=0 
    for j = 1: max(size(ml))
        M(ml(j), :) = zeros(1, params.n);
        M(:, ml(j)) = zeros(params.n ,1);
        M(ml(j), ml(j)) =1;
    end
    
end

if max(size(mr)) ~=0
    ir = params.n - 4;
    for k = 1:max(size(mr))
        M(mr(k)+ir, :) = zeros(1, params.n);
        M(:, mr(k)+ir) = zeros(params.n ,1);
        M(mr(k)+ir, mr(k)+ir) = 1;
    end
end 


end