A = 2*muc*KY2*(2*mux-Czadot);
B = -2*muc*KY2*Cza - (2*muc+Czq)*Cmadot - (2*muc-Czadot)*Cmq;
C = Cza*Cmq - (2*muc+Czq)*CMa;
z = 1+2i
[eig_val_1, eig_val_2] = (-B+sqrt(4*A*C-B^2)j)/(2*A), (-B-sqrt(4*A*C-B^2)j)/(2*A)