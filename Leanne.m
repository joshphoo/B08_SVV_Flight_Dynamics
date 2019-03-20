
figure(2);
lsim(sys_s2,de,t,x0)
title('Short Period')


[ya1,ta1,xa1] = lsim(sys_a1,dar,t,x0);
ya1(:,2) = ya1(:,2) + phi0;
ya1(:,3) = ya1(:,3) + p;
ya1(:,4) = ya1(:,4) + r;

[ya2,ta2,xa2] = lsim(sys_a2,dar,t,x0);
ya2(:,2) = ya2(:,2) + phi0;
ya2(:,3) = ya2(:,3) + p;
ya2(:,4) = ya2(:,4) + r;
[ya3,ta3,xa3] = lsim(sys_a3,dar,t,x0);
ya3(:,2) = ya3(:,2) + phi0;
ya3(:,3) = ya3(:,3) + p;
ya3(:,4) = ya3(:,4) + r;
[ya4,ta4,xa4] = lsim(sys_a4,dar,t,x0);
ya4(:,2) = ya4(:,2) + phi0;
ya4(:,3) = ya4(:,3) + p;
ya4(:,4) = ya4(:,4) + r;
    CL = 2*W/(rho*V0^2*S);
    CX0    = W*sin(th0)/(0.5*rho*V0^2*S);
    CZ0    = -W*cos(th0)/(0.5*rho*V0^2*S);
    C(1,1) = C(1,1)*V0;
    C(4,4) = C(4,4)*V0/c;
    
    C1_symmetric = [-2*muc*c/V0, 0, 0, 0; 0, (CZadot - 2*muc)*c/V0, 0, 0; 0, 0, -c/V0, 0; 0, Cmadot*c/V0, 0, -2*muc*KY2*(c/V0)];
    C2_symmetric = [CXu, CXa, CZ0, CXq; CZu, CZa, -CX0, (CZq+2*muc); 0, 0, 0, 1; Cmu Cma 0 Cmq];
    C3_symmetric = [CXde; CZde;0;Cmde];
    A1_symmetric = inv(-C1_symmetric)*C2_symmetric;
    B1_symmetric = inv(-C1_symmetric)*C3_symmetric;
    eig_symmetric = eig(A1_symmetric);
    sys_s = ss(A1_symmetric,B1_symmetric,C,D,'StateName',{'Velocity' 'Angle of attack' 'Pitch angle' 'Pitch rate'}, 'StateUnit', {'m/s' 'rad' 'rad' 'rad/s'}, 'OutputName',{'Velocity' 'Angle of attack' 'Pitch angle' 'Pitch rate'}, 'OutputUnit', {'m/s' 'rad' 'rad' 'rad/s'});
end

function [sys_a,eig_asymmetric] = ss_a(V0,hp0,m,rho0,lambda,Temp0,g,R,S,b,CYbdot,KX2,KXZ,KZ2,Cnbdot,CYb,CYp,CYr,Clb,Clp,Clr,Cnb,Cnp,Cnr,CYda,CYdr,Clda,Cldr,Cnda,Cndr,C,D)
    rho = rho0*((1+(lambda*hp0/Temp0)))^(-((g/(lambda*R))+1));   
    W = m*g;                                                
    mub = m/(rho*S*b);
    CL = 2*W/(rho*V0^2*S);
    C(3,3) = C(3,3)*2*V0/b;
    C(4,4) = C(4,4)*2*V0/b;
    
    C1_asymmetric = [(CYbdot-2*mub)*b/V0, 0, 0 ,0; 0, -b/(2*V0), 0, 0; 0, 0, -4*mub*KX2*b/V0, 4*mub*KXZ*b/V0; Cnbdot*b/V0, 0, 4*mub*KXZ*b/V0, -4*mub*KZ2*b/V0];
    C2_asymmetric = [CYb, CL, CYp, (CYr-4*mub);0, 0, 1, 0; Clb, 0, Clp, Clr; Cnb, 0 , Cnp, Cnr];
    C3_asymmetric = [CYda,CYdr; 0, 0; Clda, Cldr; Cnda, Cndr];
    A1_asymmetric = -inv(C1_asymmetric)*C2_asymmetric;
    B1_asymmetric = -inv(C1_asymmetric)*C3_asymmetric;
    eig_asymmetric = eig(A1_asymmetric);
    sys_a = ss(A1_asymmetric,B1_asymmetric,C,D,'StateName',{'Sideslip angle' 'Roll angle' 'Roll rate' 'Yaw rate'}, 'StateUnit', {'rad' 'rad' 'rad/s' 'rad/s'}, 'OutputName',{'Sideslip angle' 'Roll angle' 'Roll rate' 'Yaw rate'}, 'OutputUnit', {'rad' 'rad' 'rad/s' 'rad/s'});
end

