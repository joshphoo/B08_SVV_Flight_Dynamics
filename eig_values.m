A = 2*muc*KY2*(2*mux-Czadot);
B = -2*muc*KY2*Cza - (2*muc+Czq)*Cmadot - (2*muc-Czadot)*Cmq;
C = Cza*Cmq - (2*muc+Czq)*CMa;

[eig_val_1, eig_val_2] = (-B+sqrt(4*A*C-B^2)j)/(2*A), (-B-sqrt(4*A*C-B^2)j)/(2*A)


A_eig = 2*muc*KY2*(2*muc-CZadot)
B_eig = -2*muc*KY2*CZa - (2*muc+CZq)*Cmadot - (2*muc-CZadot)*Cmq
C_eig = CZa*Cmq - (2*muc+CZq)*Cma
eig_value_1 = (-B_eig+sqrt(4*A_eig*C_eig-B_eig^2)*i)/(2*A_eig)
eig_value_2 = (-B_eig-sqrt(4*A_eig*C_eig-B_eig^2)*i)/(2*A_eig)