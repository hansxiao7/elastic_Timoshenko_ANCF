function C = elastic_c_matrix(params)

E = params.E;
v = params.v;

C = E/(1-v^2)*[1 v 0;v 1 0;0 0 (1-v)/2];

end