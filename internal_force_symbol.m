function Q_int = internal_force_symbol(e, L, params)

Q_int = zeros(12,1);

E = params.E;
v = params.v;

e1 = e(1);
e2 = e(2);
e3 = e(3);
e4 = e(4);
e5 = e(5);
e6 = e(6);
e7 = e(7);
e8 = e(8);
e9 = e(9);
e10 = e(10);
e11 = e(11);
e12 = e(12);

h = params.h;

[x5, w5] = gauss_points(5);
y5 = x5;

for i = 1:5
    x = (x5(i)+1)*L/2;
    for j = 1:5
        y = y5(j)*h/2;
        Q = Qint(E,L,e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,v,x,y);
        Q_int = Q_int + params.b * w5(i)*w5(j)* L/2 * h/2 * Q;
    end
end

end
        
        