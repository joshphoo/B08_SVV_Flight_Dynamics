%close all;
%clear;

% Constant values concerning aircraft inertia
KX2    = 0.019;
KZ2    = 0.042;
KXZ    = 0.002;
KY2    = 1.25*1.114;
 
% Aerodynamic constants
Cmac   = 0;                     % Moment coefficient about the aerodynamic centre [ ]
CNwa   = CLa;   		        % Wing normal force slope [ ]
CNha   = 2*pi*Ah/(Ah+2);        % Stabiliser normal force slope [ ]
depsda = 4/(A+2);               % Downwash gradient [ ]
 
% Stabiblity derivatives
CXu    = -0.1; %aangepast
CXa    = 0.0; %aangepast
CZu    = -0.6; %aangepast
CZa    = -6; %aangepast
Cmu    = +0.060; %aangepast

CXadot = +0.08330;
CXq    = -0.28170;
CXde   = -0.03728;


CZadot = -0.00350;
CZq    = -5.66290;
CZde   = -0.69612;


Cmadot = +0.17800;
Cmq    = -8.79415;

Cnr    =  -0.25;
Clp    = -0.55;
Cnb    =  +0.14;
Clr    = +0.21;

CYb    = -0.7500;
CYbdot =  0     ;
CYp    = -0.0304;
CYr    = +0.8495;
CYda   = -0.0400;
CYdr   = +0.2300;

Clb    = -0.10260;


Clda   = -0.23088;
Cldr   = +0.03440;


Cnbdot =   0     ;
Cnp    =  -0.0602;

Cnda   =  -0.0120;
Cndr   =  -0.0939;


%% New!!!!!!!!!!

m_max = (9165+4050)*0.453592+95+92+74+66+61+75+78+86+68;
%m_max = (9165+2750)*0.453592+92+89+84+47+68+69+73+78+91;


C = eye(4);
D = 0;


% %% Stability
% hp0 = 1000;
% V0 = 100;
% m = m_max;
% th0 = 0;
% 
% t = 0:0.1:500;
% de = 0.0*ones(length(t),1);
% 
% x0 = [0,10/180*pi(),0,0];
% 
% symmetric = @ss_s;
% [sys_t1] = symmetric(V0,hp0,m,rho0,lambda,Temp0,g,R,S,c,CZadot,Cmadot,KY2,CXu,CXa,CXq,CZu,CZa,CZq,Cmu,Cma,Cmq,CXde,CZde,Cmde,C,D,th0);
% [yt1,tt1,xt1] = lsim(sys_t1,de,t,x0);
% 
% yt1(:,1) = yt1(:,1) + V0;
% 
% % Assym
% hp0 = 1000;
% V0 = 100;
% m = m_max;
% 
% t = 0:0.1:100;
% de = 0.0*ones(length(t),2);
% 
% x0 = [10/180*pi(),0/180*pi(),0,0];
% x1 = [0/180*pi(),10/180*pi(),0,0];
% x2 = [10/180*pi(),10/180*pi(),0,0];
% asymmetric = @ss_a;
% sys_t2 = ss_a(V0,hp0,m,rho0,lambda,Temp0,g,R,S,b,CYbdot,KX2,KXZ,KZ2,Cnbdot,CYb,CYp,CYr,Clb,Clp,Clr,Cnb,Cnp,Cnr,CYda,CYdr,Clda,Cldr,Cnda,Cndr,C,D);
% [yt20,tt20,xt20] = lsim(sys_t2,de,t,x0);
% [yt21,tt21,xt21] = lsim(sys_t2,de,t,x1);
% [yt22,tt22,xt22] = lsim(sys_t2,de,t,x2);
% 
% 
% 
% figure(1);
% subplot(2,2,1);
% hold on;
% plot(tt1,yt1(:,1));
% ylabel('TAS (m/s)')
% xlabel('Time (s)')
% subplot(2,2,2);
% hold on;
% plot(tt1,yt1(:,2));
% ylabel('AOA (rad)')
% xlabel('Time (s)')
% subplot(2,2,3)
% hold on;
% plot(tt1,yt1(:,3));
% ylabel('Pitch (rad)')
% xlabel('Time (s)')
% subplot(2,2,4)
% hold on;
% plot(tt1,yt1(:,4));
% ylabel('Pitch rate (rad/s)')
% xlabel('Time (s)')
% 
% figure(2);
% subplot(2,2,1);
% hold on;
% plot(tt1(1:100),yt1(1:100,1));
% ylabel('TAS (m/s)')
% xlabel('Time (s)')
% subplot(2,2,2);
% hold on;
% plot(tt1(1:100),yt1(1:100,2));
% ylabel('AOA (rad)')
% xlabel('Time (s)')
% subplot(2,2,3)
% hold on;
% plot(tt1(1:100),yt1(1:100,3));
% ylabel('Pitch (rad)')
% xlabel('Time (s)')
% subplot(2,2,4)
% hold on;
% plot(tt1(1:100),yt1(1:100,4));
% ylabel('Pitch rate (rad/s)')
% xlabel('Time (s)')
% 
% figure(3);
% subplot(2,2,1);
% hold on;
% plot(tt20,yt20(:,1));
% ylabel('Sideslip (rad)')
% xlabel('Time (s)')
% subplot(2,2,2);
% hold on;
% plot(tt20,yt20(:,2));
% ylabel('Roll (rad)')
% xlabel('Time (s)')
% subplot(2,2,3)
% hold on;
% plot(tt20,yt20(:,3));
% ylabel('Roll Rate (rad)')
% xlabel('Time (s)')
% subplot(2,2,4)
% hold on;
% plot(tt20,yt20(:,4));
% ylabel('Yaw Rate (rad/s)')
% xlabel('Time (s)')
% 
% figure(4);
% subplot(2,2,1);
% hold on;
% plot(tt21,yt21(:,1));
% ylabel('Sideslip (rad)')
% xlabel('Time (s)')
% subplot(2,2,2);
% hold on;
% plot(tt21,yt21(:,2));
% ylabel('Roll (rad)')
% xlabel('Time (s)')
% subplot(2,2,3)
% hold on;
% plot(tt21,yt21(:,3));
% ylabel('Roll Rate (rad)')
% xlabel('Time (s)')
% subplot(2,2,4)
% hold on;
% plot(tt21,yt21(:,4));
% ylabel('Yaw Rate (rad/s)')
% xlabel('Time (s)')
% 
% figure(5);
% subplot(2,2,1);
% hold on;
% plot(tt22,yt22(:,1));
% ylabel('Sideslip (rad)')
% xlabel('Time (s)')
% subplot(2,2,2);
% hold on;
% plot(tt22,yt22(:,2));
% ylabel('Roll (rad)')
% xlabel('Time (s)')
% subplot(2,2,3)
% hold on;
% plot(tt22,yt22(:,3));
% ylabel('Roll Rate (rad)')
% xlabel('Time (s)')
% subplot(2,2,4)
% hold on;
% plot(tt22,yt22(:,4));
% ylabel('Yaw Rate (rad/s)')
% xlabel('Time (s)')

%% Phugoid
% get all of these from flight data measured
hp0    = flightdata.Dadc1_alt.data(motion.idx1)*0.3048;      	  % pressure altitude in the stationary flight condition [m]
V0     = flightdata.Dadc1_tas.data(motion.idx1)*0.514444;          % true airspeed in the stationary flight condition [m/sec]
alpha0 = flightdata.vane_AOA.data(motion.idx1)/180*pi();    % angle of attack in the stationary flight condition [rad]
th0    = flightdata.Ahrs1_Pitch.data(motion.idx1)/180*pi();       	  % pitch angle in the stationary flight condition [rad]
m      = m_max-(flightdata.rh_engine_FU.data(motion.idx1)+flightdata.lh_engine_FU.data(motion.idx1))*0.453592;           % mass [kg]
t = flightdata.time.data(motion.idx1:motion.idxe1);
t = t-t(1);
de = (flightdata.delta_e.data(motion.idx1:motion.idxe1)-(flightdata.delta_e.data(motion.idx1)))/180*pi();

symmetric = @ss_s;
[sys_s1,eig_s1,muc_s1,CZ0] = symmetric(V0,hp0,m,rho0,lambda,Temp0,g,R,S,c,CZadot,Cmadot,KY2,CXu,CXa,CXq,CZu,CZa,CZq,Cmu,Cma,Cmq,CXde,CZde,Cmde,C,D,th0);

[ys1,ts1,xs1] = lsim(sys_s1,de,t);
ys1(:,1) = ys1(:,1) + V0;
ys1(:,2) = ys1(:,2) + alpha0;
ys1(:,3) = ys1(:,3) + th0;

A_eig_s1 = 2*muc_s1*(CZa*Cmq-2*muc_s1*Cma);
B_eig_s1 = 2*muc_s1*(CXu*Cma-Cmu*CXa)+Cmq*(CZu*CXa-CXu*CZa);
C_eig_s1 = CZ0*(Cmu*CZa-CZu*Cma);
eig_value_1_s1 = (-B_eig_s1+sqrt(4*A_eig_s1*C_eig_s1-B_eig_s1^2)*1i)/(2*A_eig_s1)*V0/c;
eig_value_2_s1 = (-B_eig_s1-sqrt(4*A_eig_s1*C_eig_s1-B_eig_s1^2)*1i)/(2*A_eig_s1)*V0/c;

%% Short Period
% input
hp0    = flightdata.Dadc1_alt.data(motion.idx2)*0.3048;      	  % pressure altitude in the stationary flight condition [m]
V0     = flightdata.Dadc1_tas.data(motion.idx2)*0.514444;          % true airspeed in the stationary flight condition [m/sec]
alpha0 = flightdata.vane_AOA.data(motion.idx2)/180*pi();    % angle of attack in the stationary flight condition [rad]
th0    = flightdata.Ahrs1_Pitch.data(motion.idx2)/180*pi();       	  % pitch angle in the stationary flight condition [rad]
m      = m_max-(flightdata.rh_engine_FU.data(motion.idx2)+flightdata.lh_engine_FU.data(motion.idx2))*0.453592;           % mass [kg]

t = flightdata.time.data(motion.idx2:motion.idxe2);
t = t-t(1);
de = (flightdata.delta_e.data(motion.idx2:motion.idxe2)-(flightdata.delta_e.data(motion.idx2)))/180*pi();

symmetric = @ss_s;
[sys_s2,eig_s2,muc_s2,CZ0] = symmetric(V0,hp0,m,rho0,lambda,Temp0,g,R,S,c,CZadot,Cmadot,KY2,CXu,CXa,CXq,CZu,CZa,CZq,Cmu,Cma,Cmq,CXde,CZde,Cmde,C,D,th0);

[ys2,ts2,xs2] = lsim(sys_s2,de,t);
ys2(:,1) = ys2(:,1) + V0;
ys2(:,2) = ys2(:,2) + alpha0;
ys2(:,3) = ys2(:,3) + th0;

A_eig = 2*muc_s2*KY2*(2*muc_s2-CZadot);
B_eig = -2*muc_s2*KY2*CZa - (2*muc_s2+CZq)*Cmadot - (2*muc_s2-CZadot)*Cmq;
C_eig = CZa*Cmq - (2*muc_s2+CZq)*Cma;
eig_value_1_s2 = (-B_eig+sqrt(4*A_eig*C_eig-B_eig^2)*1i)/(2*A_eig)*V0/c;
eig_value_2_s2 = (-B_eig-sqrt(4*A_eig*C_eig-B_eig^2)*1i)/(2*A_eig)*V0/c;

%% Dutch Roll
% input
hp0    = flightdata.Dadc1_alt.data(motion.idx3)*0.3048;      	  % pressure altitude in the stationary flight condition [m]
V0     = flightdata.Dadc1_tas.data(motion.idx3)*0.514444;          % true airspeed in the stationary flight condition [m/sec]
m      = m_max-(flightdata.rh_engine_FU.data(motion.idx3)+flightdata.lh_engine_FU.data(motion.idx3))*0.453592;           % mass [kg]
phi0   = flightdata.Ahrs1_Roll.data(motion.idx3)/180*pi();
p      = (pi()/180)*flightdata.Ahrs1_bRollRate.data(motion.idx3);
r      = (pi()/180)*flightdata.Ahrs1_bYawRate.data(motion.idx3);

t = flightdata.time.data(motion.idx3:motion.idxe3);
t = t-t(1);
dar = [];
dar(:,1) = (flightdata.delta_a.data(motion.idx3:motion.idxe3)-(flightdata.delta_a.data(motion.idx3)))/180*pi();
dar(:,2) = -(flightdata.delta_r.data(motion.idx3:motion.idxe3)-(flightdata.delta_r.data(motion.idx3)))/180*pi();

% calc
asymmetric = @ss_a;
[sys_a1,eig_a1,mub] = ss_a(V0,hp0,m,rho0,lambda,Temp0,g,R,S,b,CYbdot,KX2,KXZ,KZ2,Cnbdot,CYb,CYp,CYr,Clb,Clp,Clr,Cnb,Cnp,Cnr,CYda,CYdr,Clda,Cldr,Cnda,Cndr,C,D);

[ya1,ta1,xa1] = lsim(sys_a1,dar,t);
ya1(:,2) = ya1(:,2) + phi0;
ya1(:,3) = ya1(:,3) + p;
ya1(:,4) = ya1(:,4) + r;

eig_value_a1_1 = (2*(Cnr+2*KZ2*CYb)+sqrt(64*KZ2*(4*mub*Cnb+CYb*Cnr)-4*(Cnr+2*KZ2*CYb)^2)*1i)/(16*mub*KZ2)*V0/b;
eig_value_a1_2 = (2*(Cnr+2*KZ2*CYb)-sqrt(64*KZ2*(4*mub*Cnb+CYb*Cnr)-4*(Cnr+2*KZ2*CYb)^2)*1i)/(16*mub*KZ2)*V0/b;

%% Dutch Roll Damped
% input
hp0    = flightdata.Dadc1_alt.data(motion.idx4)*0.3048;      	  % pressure altitude in the stationary flight condition [m]
V0     = flightdata.Dadc1_tas.data(motion.idx4)*0.514444;          % true airspeed in the stationary flight condition [m/sec]
m      = m_max-(flightdata.rh_engine_FU.data(motion.idx4)+flightdata.lh_engine_FU.data(motion.idx4))*0.453592;           % mass [kg]
phi0   = flightdata.Ahrs1_Roll.data(motion.idx4)/180*pi();
p      = (pi/180)*flightdata.Ahrs1_bRollRate.data(motion.idx4);
r      = (pi/180)*flightdata.Ahrs1_bYawRate.data(motion.idx4);

t = flightdata.time.data(motion.idx4:motion.idxe4);
t = t-t(1);
dar = [];
dar(:,1) = (flightdata.delta_a.data(motion.idx4:motion.idxe4)-(flightdata.delta_a.data(motion.idx6)))/180*pi();
dar(:,2) = -(flightdata.delta_r.data(motion.idx4:motion.idxe4)-(flightdata.delta_r.data(motion.idx6)))/180*pi();

% calc
asymmetric = @ss_a;
[sys_a2,eig_a2,mub] = ss_a(V0,hp0,m,rho0,lambda,Temp0,g,R,S,b,CYbdot,KX2,KXZ,KZ2,Cnbdot,CYb,CYp,CYr,Clb,Clp,Clr,Cnb,Cnp,Cnr,CYda,CYdr,Clda,Cldr,Cnda,Cndr,C,D);

[ya2,ta2,xa2] = lsim(sys_a2,dar,t);
ya2(:,2) = ya2(:,2) + phi0;
ya2(:,3) = ya2(:,3) + p;
ya2(:,4) = ya2(:,4) + r;

%% Aperiodic Roll
% input
hp0    = flightdata.Dadc1_alt.data(motion.idx5)*0.3048;      	  % pressure altitude in the stationary flight condition [m]
V0     = flightdata.Dadc1_tas.data(motion.idx5)*0.514444;          % true airspeed in the stationary flight condition [m/sec]
m      = m_max-(flightdata.rh_engine_FU.data(motion.idx5)+flightdata.lh_engine_FU.data(motion.idx5))*0.453592;           % mass [kg]
phi0   = - flightdata.Ahrs1_Roll.data(motion.idx5)/180*pi();          % mass [kg]
p      = -(pi/180)*flightdata.Ahrs1_bRollRate.data(motion.idx5);
r      = -(pi/180)*flightdata.Ahrs1_bYawRate.data(motion.idx5);


t = flightdata.time.data(motion.idx5:motion.idxe5);
t = t-t(1);
dar = [];
dar(:,1) = (flightdata.delta_a.data(motion.idx5:motion.idxe5)-flightdata.delta_a.data(motion.idx5))/180*pi();
dar(:,2) = -(flightdata.delta_r.data(motion.idx5:motion.idxe5)-flightdata.delta_r.data(motion.idx5))/180*pi();

% calc
asymmetric = @ss_a;
[sys_a3,eig_a3,mub] = ss_a(V0,hp0,m,rho0,lambda,Temp0,g,R,S,b,CYbdot,KX2,KXZ,KZ2,Cnbdot,CYb,CYp,CYr,Clb,Clp,Clr,Cnb,Cnp,Cnr,CYda,CYdr,Clda,Cldr,Cnda,Cndr,C,D);

[ya3,ta3,xa3] = lsim(sys_a3,dar,t);
ya3(:,2) = ya3(:,2) + phi0;
ya3(:,3) = ya3(:,3) + p;
ya3(:,4) = ya3(:,4) + r;

eig_value_a3_1 = Clp/(4*mub*KX2)*V0/b;

%% Spiral
% input
hp0    = flightdata.Dadc1_alt.data(motion.idx6)*0.3048;      	  % pressure altitude in the stationary flight condition [m]
V0     = flightdata.Dadc1_tas.data(motion.idx6)*0.514444;          % true airspeed in the stationary flight condition [m/sec]
m      = m_max-(flightdata.rh_engine_FU.data(motion.idx6)+flightdata.lh_engine_FU.data(motion.idx6))*0.453592;           % mass [kg]
phi0   = -flightdata.Ahrs1_Roll.data(motion.idx6)/180*pi();          % mass [kg]
p      = -(pi/180)*(flightdata.Ahrs1_bRollRate.data(motion.idx6));
r      = -(pi/180)*(flightdata.Ahrs1_bYawRate.data(motion.idx6));


t = flightdata.time.data(motion.idx6:motion.idxe6);
t = t-t(1);
dar = [];
dar(1,:) = (flightdata.delta_a.data(motion.idx6:motion.idxe6)-(flightdata.delta_a.data(motion.idx6)))/180*pi();
dar(2,:) = -(flightdata.delta_r.data(motion.idx6:motion.idxe6)-(flightdata.delta_r.data(motion.idx6)))/180*pi();

% calc
asymmetric = @ss_a;
[sys_a4,eig_a4,mub,CL] = ss_a(V0,hp0,m,rho0,lambda,Temp0,g,R,S,b,CYbdot,KX2,KXZ,KZ2,Cnbdot,CYb,CYp,CYr,Clb,Clp,Clr,Cnb,Cnp,Cnr,CYda,CYdr,Clda,Cldr,Cnda,Cndr,C,D);

[ya4,ta4,xa4] = lsim(sys_a4,dar,t);
ya4(:,2) = ya4(:,2) + phi0;
ya4(:,3) = ya4(:,3) + p;
ya4(:,4) = ya4(:,4) + r;

eig_value_a4_1 = (2*CL*(Clb*Cnr-Cnb*Clr))/(Clp*(CYb*Cnr+4*mub*Cnb)-Cnp*(CYb*Clr+4*mub*Clb))*V0/b;

% Plotten
figure(1)
subplot(2,2,1);
hold on;
plot(ts1,ys1(:,1));
legend({'measured','calculated'},'Location','northeast');
legend('boxoff');

subplot(2,2,2);
hold on;
plot(ts1,ys1(:,2));
legend({'measured','calculated'},'Location','northeast');
legend('boxoff');
subplot(2,2,3)
hold on;
plot(ts1,ys1(:,3));
legend({'measured','calculated'},'Location','northeast');
legend('boxoff');
subplot(2,2,4)
hold on;
plot(ts1,ys1(:,4));
legend({'measured','calculated'},'Location','southeast');
legend('boxoff');

figure(2)
subplot(2,2,1);
hold on;
plot(ts2,ys2(:,1));
legend({'measured','calculated'},'Location','northeast');
legend('boxoff');
subplot(2,2,2);
hold on;
plot(ts2,ys2(:,2));
legend({'measured','calculated'},'Location','east');
legend('boxoff');
subplot(2,2,3)
hold on;
plot(ts2,ys2(:,3));
legend({'measured','calculated'},'Location','south');
legend('boxoff');
subplot(2,2,4)
hold on;
plot(ts2,ys2(:,4));
legend({'measured','calculated'},'Location','northeast');
legend('boxoff');

figure(3)
subplot(2,2,1);
hold on;
plot(ta1,ya1(:,2));
legend({'measured','calculated'},'Location','southeast');
legend('boxoff');
subplot(2,2,2);
hold on;
plot(ta1,ya1(:,3));
legend({'measured','calculated'},'Location','southeast');
legend('boxoff');
subplot(2,2,3)
hold on;
plot(ta1,ya1(:,4));
legend({'measured','calculated'},'Location','southeast');
legend('boxoff');

figure(4)
subplot(2,2,1);
hold on;
plot(ta2,ya2(:,2));
legend({'measured','calculated'},'Location','southeast');
legend('boxoff');
subplot(2,2,2);
hold on;
plot(ta2,ya2(:,3));
legend({'measured','calculated'},'Location','southeast');
legend('boxoff');
subplot(2,2,3)
hold on;
plot(ta2,ya2(:,4));
legend({'measured','calculated'},'Location','southeast');
legend('boxoff');

figure(5)
subplot(2,2,1);
hold on;
plot(ta3,ya3(:,2));
legend({'measured','calculated'},'Location','southeast');
legend('boxoff');
subplot(2,2,2);
hold on;
plot(ta3,ya3(:,3));
legend({'measured','calculated'},'Location','southwest');
legend('boxoff');
subplot(2,2,3)
hold on;
plot(ta3,ya3(:,4));
legend({'measured','calculated'},'Location','southeast');
legend('boxoff');

figure(6)
subplot(2,2,1);
hold on;
plot(ta4,ya4(:,2));
legend({'measured','calculated'},'Location','northwest');
legend('boxoff');
subplot(2,2,2);
hold on;
plot(ta4,ya4(:,3));
legend({'measured','calculated'},'Location','southeast');
legend('boxoff');
subplot(2,2,3)
hold on;
plot(ta4,ya4(:,4));
legend({'measured','calculated'},'Location','northwest');
legend('boxoff');
hold off;

savefig(refmotion,'realmotions.fig')
saveas(refmotion(1),'realmotions1Phugoid.png')
saveas(refmotion(2),'realmotions2ShortPeriod.png')
saveas(refmotion(3),'realmotions3Droll.png')
saveas(refmotion(4),'realmotions4YDDroll.png')
saveas(refmotion(5),'realmotions5Aroll.png')
saveas(refmotion(5),'realmotions6spiral.png')


fprintf( 'Symmetric: Analytical Eigenvalues'), disp(0);
fprintf( 'Eigenvalue Short Period'), disp(eig_value_1_s2);
fprintf( 'Eigenvalue Short Period'), disp(eig_value_2_s2);
fprintf( 'Eigenvalue Phugoid'), disp(eig_value_1_s1);
fprintf( 'Eigenvalue Phugoid'), disp(eig_value_2_s1);
fprintf( 'Symmetric: Numerical Eigenvalues'), disp(eig_s1);


fprintf( 'Asymmetric: Analytical Eigenvalues'), disp(0);
fprintf( 'Eigenvalue Aperiodic Roll'), disp(eig_value_a3_1);
fprintf( 'Eigenvalue Dutch Roll'), disp(eig_value_a1_1);
fprintf( 'Eigenvalue Dutch Roll'), disp(eig_value_a1_2);
fprintf( 'Eigenvalue Spiral'), disp(eig_value_a4_1);
fprintf( 'Asymmetric: Numerical Eigenvalues'), disp(eig_a1);


% ys1_1n = [];
% ys1_1n(:,1) = ys1(find(ys1(:,1) == min(ys1(:,1))):end,1)-mean(ys1(:,2));
% idx0 = [];
% for i = 1:length(ys1_1n(:,1))
%     if i ~= 1
%         if sign(ys1_1n(i,1))+sign(ys1_1n(i-1,1)) == 0
%             idx0 = [idx0 ; find(ys1_1n(:,1) == ys1_1n(i,1))];
%         end
%     end
% end
% figure(1)
% P = [];
% for i = 3:length(idx0)
% P(i-2) = (idx0(i)-idx0(i-2))/10;
% end
% ys1_max = 0;
% for j = 1:(length(idx0)-1)/2
%     ys1_max(j) = max(ys1_1n(idx0(2*j-1):idx0(2*j),1));
%     max_idx(j) = find(ys1_1n(:,1) == max(ys1_1n(idx0(2*j-1):idx0(2*j),1)));
% end
% ys1_maxc = ys1_max/ys1_max(1)
% 
% ys1_min = 0;
% ys1_min(1) = ys1_1n(1);
% min_idx(1) = 1;
% for j = 1:(length(idx0)-1)/2
%     ys1_min(j+1) = min(ys1_1n(idx0(2*j):idx0(2*j+1),1));
%     min_idx(j+1) = find(ys1_1n(:,1) == min(ys1_1n(idx0(2*j):idx0(2*j+1),1)))
% end
% ys1_minc = ys1_min/ys1_min(1)
% 
% 
% for i = 1:length(ys1_maxc)
% if ys1_maxc(i) < 0.5
%     half = ((0.5-ys1_maxc(i))*max_idx(i-1) + (ys1_maxc(i-1)-0.5)*max_idx(i))/(ys1_maxc(i-1)-ys1_maxc(i))/10
%     break
% end
% end
% 
% for i = 1:length(ys1_minc)
% if ys1_minc(i) < 0.5
%     half = ((0.5-ys1_minc(i))*min_idx(i-1) + (ys1_minc(i-1)-0.5)*min_idx(i))/(ys1_minc(i-1)-ys1_minc(i))/10
%     break
% end
% end
% 
% Pavg = mean(P);
% plot(1:length(ys1_1n(:,1)),ys1_1n(:,1))


% ya1_1n = [];
% ya1_1n(:,2) = ya1(:,2)-mean(ya1(:,2));
% idx0 = [];
% for i = 1:length(ya1_1n(:,1))
%     if i ~= 1
%         if sign(ya1_1n(i,1))+sign(ya1_1n(i-1,1)) == 0
%             idx0 = [idx0 ; find(ya1_1n(:,1) == ya1_1n(i,1))];
%         end
%     end
% end

% P = [];
% for i = 3:length(idx0)
% P(i-2) = (idx0(i)-idx0(i-2))/10;
% end
% ys3_max = 0;
% for j = 1:(length(idx0)-1)/2
%     ys3_max(j) = max(ys3_1n(idx0(2*j-1):idx0(2*j),1));
%     max_idx(j) = find(ys3_1n(:,1) == max(ys3_1n(idx0(2*j-1):idx0(2*j),1)));
% end
% ys3_maxc = ys3_max/ys3_max(1)
% 
% ys3_min = 0;
% ys3_min(1) = ys3_1n(1);
% min_idx(1) = 1;
% for j = 1:(length(idx0)-1)/2
%     ys3_min(j+1) = min(ys3_1n(idx0(2*j):idx0(2*j+1),1));
%     min_idx(j+1) = find(ys3_1n(:,1) == min(ys3_1n(idx0(2*j):idx0(2*j+1),1)))
% end
% ys3_minc = ys3_min/ys3_min(1)
% 
% 
% for i = 1:length(ys3_maxc)
% if ys3_maxc(i) < 0.5
%     half2 = ((0.5-ys3_maxc(i))*max_idx(i-1) + (ys3_maxc(i-1)-0.5)*max_idx(i))/(ys3_maxc(i-1)-ys3_maxc(i))/10
%     break
% end
% end
% 
% for i = 1:length(ys3_minc)
% if ys3_minc(i) < 0.5
%     half2 = ((0.5-ys3_minc(i))*min_idx(i-1) + (ys3_minc(i-1)-0.5)*min_idx(i))/(ys3_minc(i-1)-ys3_minc(i))/10
%     break
% end
% end
% 
% Pavg = mean(P);
% plot(1:length(ya1_1n(:,1)),ya1_1n(:,1))


% thalf = [];
% for h = 1:2
% vmax = [];
% vmin = [];
% vmaxi = [];
% vmini = []
% vn = [];
% vmaxc = [];
% vminc = [];
% for i = 1:length(ys1(:,h))
%     if i ~= 1 & i~= length(ys1(:,h))
%         if ys1(i,h) > ys1(i+1,h) & ys1(i,h) > ys1(i-1,h) & median(ys1(:,h) < ys1(i,h))
%         vmax = [vmax;find(ys1(:,h) == ys1(i,h))];
%         elseif ys1(i,h) < ys1(i+1,h) & ys1(i,h) < ys1(i-1,h) & median(ys1(:,h) > ys1(i,h))
%         vmin = [vmin;find(ys1(:,h) == ys1(i,h))];
%         end  
%     end
% end
% vmaxi = vmax(find(ys1(vmax,h)==max(ys1(vmax,h))):end);
% vmax = ys1(vmax(find(ys1(vmax,h)==max(ys1(vmax,h))):end),h);
% vmini = vmin(find(ys1(vmin,h)==min(ys1(vmin,h))):end);
% vmin = ys1(vmin(find(ys1(vmin,h)==min(ys1(vmin,h))):end),h);
% 
% for i = 1:length(vmaxi)
%     if i~= 1
%         if abs(vmaxi(i)-vmaxi(i-1)) < 100
%         if vmax(i) > vmax(i-1)
%             vmax(i-1) = 0;
%             vmaxi(i-1) = 0;
%         else
%             vmax(i) = 0;
%             vmaxi(i) = 0;
%             if i~= length(vmaxi)
%                 if abs(vmaxi(i+1)-vmaxi(i-1)) < 100
%                      if vmax(i+1) > vmax(i-1)
%                           %vmax(i-1) = 0;
%                           %vmaxi(i-1) = 0;
%                      else
%                           %vmax(i+1) = 0;
%                           %vmaxi(i+1) = 0;
%                      end
%                 end
%             end
%         end
%         end
%     end
% end
% 
% % vmaxi = vmaxi(vmaxi~=0);
% % vmax = vmax(vmax~=0);
% vn = ys1(vmini(1):vmaxi(end),h);
% %vmax = vmax-median(vn);
% vmin = vmin-median(vn);
% vn = vn-mean(vn);
% figure(h)
% plot(1:length(vn),vn)
% vmaxc = vmax/vmax(1);
% vminc = vmin/vmin(1);
% 
% for i = 1:length(vmaxc)
% if vmaxc(i) < 0.5
%     thalf(h,1) = ((0.5-vmaxc(i))*vmaxi(i-1)+(vmaxc(i-1)-0.5)*vmaxi(i))/(vmaxc(i-1)-vmaxc(i))/10
%     break
% end
% end
% 
% for i = 1:length(vminc)
% if vminc(i) < 0.5
%     thalf2(h) = ((0.5-vminc(i))*vmini(i-1)+(vminc(i-1)-0.5)*vmini(i))/(vminc(i-1)-vminc(i))/10
%     break
% end
% end
% end

%% Functions
function [sys_s,eig_symmetric,muc,CZ0] = ss_s(V0,hp0,m,rho0,lambda,Temp0,g,R,S,c,CZadot,Cmadot,KY2,CXu,CXa,CXq,CZu,CZa,CZq,Cmu,Cma,Cmq,CXde,CZde,Cmde,C,D,th0)
    rho = rho0*((1+(lambda*hp0/Temp0)))^(-((g/(lambda*R))+1));   
    W = m*g;                                                
    muc = m/(rho*S*c);
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

function [sys_a,eig_asymmetric,mub,CL] = ss_a(V0,hp0,m,rho0,lambda,Temp0,g,R,S,b,CYbdot,KX2,KXZ,KZ2,Cnbdot,CYb,CYp,CYr,Clb,Clp,Clr,Cnb,Cnp,Cnr,CYda,CYdr,Clda,Cldr,Cnda,Cndr,C,D)
    rho = rho0*((1+(lambda*hp0/Temp0)))^(-((g/(lambda*R))+1));   
    W = m*g;                                     
    mub = m/(rho*S*b);
    CL = 2*W/(rho*V0^2*S);
    C(3,3) = C(3,3)*2*V0/b;
    C(4,4) = C(4,4)*2*V0/b;
    
    C1_asymmetric = [(CYbdot-2*mub)*b/V0, 0, 0 ,0; 0, -b/(2*V0), 0, 0; 0, 0, -4*mub*KX2*b/V0, 4*mub*KXZ*b/V0; Cnbdot*b/V0, 0, 4*mub*KXZ*b/V0, -4*mub*KZ2*b/V0];
    C2_asymmetric = [CYb, CL, CYp, (CYr-4*mub);0, 0, 1, 0; Clb, 0, Clp, Clr; Cnb, 0 , Cnp, Cnr];
    C3_asymmetric = [CYda,CYdr; 0, 0; Clda, Cldr; Cnda, Cndr];
    A1_asymmetric = inv(-C1_asymmetric)*C2_asymmetric;
    B1_asymmetric = inv(-C1_asymmetric)*C3_asymmetric;
    eig_asymmetric = eig(A1_asymmetric);
    sys_a = ss(A1_asymmetric,B1_asymmetric,C,D,'StateName',{'Sideslip angle' 'Roll angle' 'Roll rate' 'Yaw rate'}, 'StateUnit', {'rad' 'rad' 'rad/s' 'rad/s'}, 'OutputName',{'Sideslip angle' 'Roll angle' 'Roll rate' 'Yaw rate'}, 'OutputUnit', {'rad' 'rad' 'rad/s' 'rad/s'});
end
