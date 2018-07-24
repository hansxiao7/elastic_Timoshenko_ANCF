function Q_int = total_internal_force(e, params)

Q_int = zeros(params.n, 1);

for i = 1:params.ne
    eele = e(8*i-7:8*i+4,1);
    L = params.x(i);
    Qele = internal_force(eele, L, params);
    Q_int(8*i-7:8*i+4,1) = Q_int(8*i-7:8*i+4,1) + Qele;
end

ml = params.ml;
mr = params.mr;
if max(size(ml)) ~=0 
    for j = 1: max(size(ml))
        Q_int(ml(j), 1) = 0;
    end
    
end

if max(size(mr)) ~=0
    ir = params.n - 4;
    for k = 1:max(size(mr))
        Q_int(mr(k)+ir, 1) = 0;
    end
end 


end
