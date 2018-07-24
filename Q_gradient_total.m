function Q_gradient_total = Q_gradient_total(e, params)

Q_gradient_total = zeros(params.n);
E = params.E;
v = params.v;
h = params.h;

[x5, w5] = gauss_points(5);
y5 = x5;

for i = 1: params.ne
    eele = e(8*i-7:8*i+4,1);
    L = params.x(i);
    e1 = eele(1);
    e2 = eele(2);
    e3 = eele(3);
    e4 = eele(4);
    e5 = eele(5);
    e6 = eele(6);
    e7 = eele(7);
    e8 = eele(8);
    e9 = eele(9);
    e10 = eele(10);
    e11 = eele(11);
    e12 = eele(12);
    Qe_gradient = zeros(12);
    
    for j = 1:5
        x = (x5(j)+1) * L/2;
        for k = 1:5
            y = y5(k) * h/2;
            Qe_gradient = Qe_gradient + params.b * w5(j)*w5(k)* L/2 * h/2 * Q_gradient(E,L,e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,v,x,y);
        end
    end
    
    ii = (i-1)*8 + 1;
    
    Q_gradient_total(ii:ii+11, ii:ii+11) = Q_gradient_total(ii:ii+11, ii:ii+11) + Qe_gradient;
    
end

end

            