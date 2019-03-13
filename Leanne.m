
% Citation 550 - Linear simulation

% xcg = 0.25*c

% Stationary flight condition

hp0    = 1500;      	  % pressure altitude in the stationary flight condition [m]
V0     = 75;          % true airspeed in the stationary flight condition [m/sec]
alpha0 = 5/180*pi;    % angle of attack in the stationary flight condition [rad]
th0    = 0;       	  % pitch angle in the stationary flight condition [rad]

% Aircraft mass
m      = 6000;           % mass [kg]

% aerodynamic properties
e      = 0.8;           % Oswald factor [ ]
CD0    = 0.04;          % Zero lift drag coefficient [ ]
CLa    = 5.084;         % Slope of CL-alpha curve [ ]

% Longitudinal stability
Cma    = -0.5626;           % longitudinal stabilty [ ]
Cmde   = -1.1642;           % elevator effectiveness [ ]

% Aircraft geometry

S      = 30.00;	          % wing area [m^2]
Sh     = 0.2*S;           % stabiliser area [m^2]
Sh_S   = Sh/S;	          % [ ]
lh     = 0.71*5.968;      % tail length [m]
c      = 2.0569;          % mean aerodynamic cord [m]
lh_c   = lh/c;	          % [ ]
b      = 15.911;          % wing span [m]
bh     = 5.791;	          % stabilser span [m]
A      = b^2/S;           % wing aspect ratio [ ]
Ah     = bh^2/Sh;         % stabiliser aspect ratio [ ]
Vh_V   = 1;               % [ ]
ih     = -2*pi/180;       % stabiliser angle of incidence [rad]

% Constant values concerning atmosphere and gravity

rho0   = 1.2250;          % air density at sea level [kg/m^3] 
lambda = -0.0065;         % temperature gradient in ISA [K/m]
Temp0  = 288.15;          % temperature at sea level in ISA [K]
R      = 287.05;          % specific gas constant [m^2/sec^2K]
g      = 9.81;            % [m/sec^2] (gravity constant)
gamma  = 1.4;             % air ratio specific heats 
p0     = 101325;          % sea level pressure [Pa]

rho    = rho0*((1+(lambda*hp0/Temp0)))^(-((g/(lambda*R))+1));   % [kg/m^3]  (air density)
W      = m*g;                                                   % [N]       (aircraft weight)

% Constant values concerning aircraft inertia

muc    = m/(rho*S*c);
mub    = m/(rho*S*b);
KX2    = 0.019;
KZ2    = 0.042;
KXZ    = 0.002;
KY2    = 1.25*1.114;

% Aerodynamic constants

Cmac   = 0;                     % Moment coefficient about the aerodynamic centre [ ]
CNwa   = CLa;   		        % Wing normal force slope [ ]
CNha   = 2*pi*Ah/(Ah+2);        % Stabiliser normal force slope [ ]
depsda = 4/(A+2);               % Downwash gradient [ ]

% Lift and drag coefficient

CL = 2*W/(rho*V0^2*S);               % Lift coefficient [ ]
CD = CD0 + (CLa*alpha0)^2/(pi*A*e);  % Drag coefficient [ ]

% Stabiblity derivatives

CX0    = W*sin(th0)/(0.5*rho*V0^2*S);
CXu    = -0.02792;
CXa    = -0.47966;
CXadot = +0.08330;
CXq    = -0.28170;
CXde   = -0.03728;

CZ0    = -W*cos(th0)/(0.5*rho*V0^2*S);
CZu    = -0.37616;
CZa    = -5.74340;
CZadot = -0.00350;
CZq    = -5.66290;
CZde   = -0.69612;

Cmu    = +0.06990;
Cmadot = +0.17800;
Cmq    = -8.79415;

CYb    = -0.7500;
CYbdot =  0     ;
CYp    = -0.0304;
CYr    = +0.8495;
CYda   = -0.0400;
CYdr   = +0.2300;

Clb    = -0.10260;
Clp    = -0.71085;
Clr    = +0.23760;
Clda   = -0.23088;
Cldr   = +0.03440;

Cnb    =  +0.1348;
Cnbdot =   0     ;
Cnp    =  -0.0602;
Cnr    =  -0.2061;
Cnda   =  -0.0120;
Cndr   =  -0.0939;

syms u alpha theta q udot alphadot thetadot qdot de
C = eye(4);
D = 0;

%symmetric
% state_vector_symmetric = [u; alpha; theta; q];
% state_vector_dot_symmetric = [udot; alphadot; thetadot; qdot];
% input_vector_symmetric = [de];

C1_symmetric = [-2*muc*c/V0^2, 0, 0, 0; 0, (CZadot - 2*muc)*c/V0, 0, 0; 0, 0, -c/V0, 0; 0, Cmadot*c/V0, 0, -2*muc*KY2*(c/V0)^2];
C2_symmetric = [CXu/V0, CXa, CZ0, CXq*c/V0; CZu/V0, CZa, -CX0, (CZq+2*muc)*c/V0; 0, 0, 0, c/V0; Cmu/V0 Cma 0 Cmq*c/V0];
C3_symmetric = [CXde; CZde;0;Cmde];
A1_symmetric = -inv(C1_symmetric)*C2_symmetric;
B1_symmetric = -inv(C1_symmetric)*C3_symmetric;
A2 = [V0*CXu/(2*muc*c), V0^2*CXa/(2*muc*c), V0^2*CZ0/(2*muc*c), V0*CXq/(2*muc); -CZu/((CZadot-2*muc)*c), -V0*CZa/((CZadot-2*muc)*c), V0*CX0/((CZadot-2*muc)*c), -(CZq+2*muc)/(CZadot-2*muc); 0, 0, 0, 1; V0/(2*muc*KY2*c^2)*(Cmu-Cmadot*CZu/(CZadot-2*muc)),V0^2/(2*muc*KY2*c^2)*(Cma-Cmadot*CZa/(CZadot-2*muc)),V0^2/(2*muc*KY2*c^2)*(Cmadot*CX0/(CZadot-2*muc)),V0/(2*muc*KY2*c)*(Cmq-Cmadot*(CZq+2*muc)/(CZadot-2*muc))];
B2 = [V0^2*CXde/(2*muc*c);-V0^2*CZde/((CZadot-2*muc)*c);0;V0^2/(2*muc*KY2*c^2)*(Cmde-(Cmadot*CZde)/(CZadot-2*muc))];

%asymmetric
C1_asymmetric = [(CYbdot-2*mub)*b/V0, 0, 0 ,0; 0, -b/(2*V0), 0, 0; 0, 0, -2*mub*KX2*b^2/V0^2, 2*mub*KXZ*b^2/V0^2; Cnb*b/V0, 0, 2*mub*KXZ*b^2/V0^2, -2*mub*KZ2*b^2/V0^2];
C2_asymmetric = [CYb, CL, CYp*b/(2*V0), (CYr-4*mub)*b/(2*V0);0, 0, b/(2*V0), 0; Clb, 0, Clp*b/(2*V0), Clr*b/(2*V0); Cnb, 0 , Cnp*b/(2*V0), Cnr*b/(2*V0)];
C3_asymmetric = [CYda,CYdr; 0, 0; Clda, Cldr; Cnda, Cndr];
A1_asymmetric = -inv(C1_asymmetric)*C2_asymmetric;
B1_asymmetric = -inv(C1_asymmetric)*C3_asymmetric;

T = 5;

sys_s = ss(A1_symmetric,B1_symmetric,C,D,'StateName',{'Velocity' 'Angle of attack' 'Pitch angle' 'Pitch rate'}, 'StateUnit', {'m/s' 'rad' 'rad' 'rad/s'}, 'OutputName',{'Velocity' 'Angle of attack' 'Pitch angle' 'Pitch rate'}, 'OutputUnit', {'m/s' 'rad' 'rad' 'rad/s'});
sys_a = ss(A1_asymmetric,B1_asymmetric,C,D);
size(sys_s)

t = 0:.1:10;
u = [0; zeros(length(t)-1,1)];
x0 = [V0; 0; th0; 0]
lsim(sys_s,u,t,x0)
isstable(sys_a)
%impulse(sys_s, 10)

% Mean0 = [V0; alpha0; th0; 0] 
% sys_s2 = ssm(A1_symmetric,B1_symmetric,C,D, 'Mean0', Mean0)
% 
% % impulse(sys_a);
% 
% y = simulate(sys_s2,T);
% fMdl = filter(sys_s2,y);
% % fDSSMMdl = filter(DSSMMdl,y);
% sMdl = smooth(sys_s2,y);
% %sDSSMMdl = smooth(DSSMMdl,y);
% % [Y,X,U,E] = simulate(sys_s2,4)
% figure;
% %plot(1:T,y,'-o',1:T,fMdl,'-d',1:T,sMdl,'-*');
% plot(1:T,y,'-o')
% %legend('Simulated Data','Filtered States -- Mdl','Smooth States -- sMdl');

% symmetric_eq1 = [CXu-2*muc*c/V0*(ddt), CXa, CZ0, CXq; CZu, CZa+(CZadot-2*muc)*c/V0*(ddt), -CX0, CZq+2*muc; 0, 0, -c/V0*(ddt), 1; Cmu, Cma+Cmadot*c/V0*(ddt),0,Cmq-2*muc*KY2]
% symmetric_eq2 = [ucirc;alpha;theta;q*c/V0]
% sym1 = symmetric_eq1*symmetric_eq2