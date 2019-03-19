
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

%rho    = rho0*((1+(lambda*hp0/Temp0)))^(-((g/(lambda*R))+1));   % [kg/m^3]  (air density)
%W      = m*g;                                                   % [N]       (aircraft weight)

% Constant values concerning aircraft inertia
%muc    = m/(rho*S*c);
%mub    = m/(rho*S*b);
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
%CL = 2*W/(rho*V0^2*S);               % Lift coefficient [ ]
%CD = CD0 + (CLa*alpha0)^2/(pi*A*e);  % Drag coefficient [ ]

% Stabiblity derivatives
%CX0    = W*sin(th0)/(0.5*rho*V0^2*S);
CXu    = -0.02792;
CXa    = -0.47966;
CXadot = +0.08330;
CXq    = -0.28170;
CXde   = -0.03728;

%CZ0    = -W*cos(th0)/(0.5*rho*V0^2*S);
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

C = eye(4);
D = 0;


%% Short period
% get all of these from flight data measured
hp0    = 250;      	  % pressure altitude in the stationary flight condition [m]
V0     = 67;          % true airspeed in the stationary flight condition [m/sec]
alpha0 = 0.1;    % angle of attack in the stationary flight condition [rad]
th0    = 0/180*pi;       	  % pitch angle in the stationary flight condition [rad]
m      = 6000;           % mass [kg]
de_input = 4/180*pi;
de_time = 1;
t = 10;
dt = 0.01;

symmetric = @ss_s;
[sys_s1,eig_s1] = symmetric(V0,hp0,m,rho0,lambda,Temp0,g,R,S,c,CZadot,Cmadot,KY2,CXu,CXa,CXq,CZu,CZa,CZq,Cmu,Cma,Cmq,CXde,CZde,Cmde,C,D,th0);

t1 = 0:dt:t;
u1 = de_input(ones(length(t1),1));
%u1 = de_input(ones(de_time/dt,1));
%u1 = [u1; zeros(t/dt-length(u1)+1,1)];
x0 = [1; alpha0; th0; 0];
figure(1);
lsim(sys_s1,u1,t1,x0)

%% Phugoid
% input
hp0    = 2500;      	  % pressure altitude in the stationary flight condition [m]
V0     = 97;          % true airspeed in the stationary flight condition [m/sec]
alpha0 = 0.1;    % angle of attack in the stationary flight condition [rad]
th0    = 0.0;       	  % pitch angle in the stationary flight condition [rad]
m      = 6000;           % mass [kg]
de_input = 0/180*pi;
de_time = 200;
t = 200;

% calc
symmetric = @ss_s;
[sys_s2,eig_s2] = symmetric(V0,hp0,m,rho0,lambda,Temp0,g,R,S,c,CZadot,Cmadot,KY2,CXu,CXa,CXq,CZu,CZa,CZq,Cmu,Cma,Cmq,CXde,CZde,Cmde,C,D,th0);

det2 = [1:de_time/dt];
t2 = 0:dt:t;
u2 = de_input(ones(length(t2),1));
u2(det2) = de_input;
% u2 = de_input(ones(de_time/dt,1));
% u2 = [u2; zeros(t/dt-length(u2)+1,1)];
x0 = [1; alpha0; th0; 0];
figure(2);
lsim(sys_s2,u2,t2,x0)



%% Aperiodic Roll
% input
hp0    = 2500;      	  % pressure altitude in the stationary flight condition [m]
V0     = 97;          % true airspeed in the stationary flight condition [m/sec]
beta0 = 0.1;
phi0 = 0.1;
th0    = 0.0;       	  % pitch angle in the stationary flight condition [rad]
m      = 6000;           % mass [kg]
da_input = 0.1/180*pi;
dr_input = 0.1/180*pi;
d_time = 200;
t = 200;

% calc
asymmetric = @ss_a;
[sys_a1,eig_a1] = ss_a(V0,hp0,m,rho0,lambda,Temp0,g,R,S,b,CYbdot,KX2,KXZ,KZ2,Cnbdot,CYb,CYp,CYr,Clb,Clp,Clr,Cnb,Cnp,Cnr,CYda,CYdr,Clda,Cldr,Cnda,Cndr,C,D,th0);

dt3 = [1:d_time/dt];
t3 = 0:dt:t;
u3 = da_input(ones(length(t3),2));
u3(dt3,1) = da_input;
u3(dt3,2) = dr_input;
% u2 = de_input(ones(de_time/dt,1));
% u2 = [u2; zeros(t/dt-length(u2)+1,1)];
x0 = [beta0; phi0; 0; 0];
figure(3);
lsim(sys_a1,u3,t3,x0)
%% Spiral
% input
hp0    = 2500;      	  % pressure altitude in the stationary flight condition [m]
V0     = 97;          % true airspeed in the stationary flight condition [m/sec]
beta0 = 0.1;
phi0 = 0.1;
th0    = 0.0;       	  % pitch angle in the stationary flight condition [rad]
m      = 6000;           % mass [kg]
da_input = 0.1/180*pi;
dr_input = 0.1/180*pi;
d_time = 200;
t = 200;

% calc
asymmetric = @ss_a;
[sys_a2,eig_a2] = ss_a(V0,hp0,m,rho0,lambda,Temp0,g,R,S,b,CYbdot,KX2,KXZ,KZ2,Cnbdot,CYb,CYp,CYr,Clb,Clp,Clr,Cnb,Cnp,Cnr,CYda,CYdr,Clda,Cldr,Cnda,Cndr,C,D,th0);

dt4 = [1:d_time/dt];
t4 = 0:dt:t;
u4 = da_input(ones(length(t4),2));
u4(dt4,1) = da_input;
u4(dt4,2) = dr_input;
% u2 = de_input(ones(de_time/dt,1));
% u2 = [u2; zeros(t/dt-length(u2)+1,1)];
x0 = [beta0; phi0; 0; 0];
figure(4);
lsim(sys_a2,u4,t4,x0)

%% Dutch Roll
% input
hp0    = 2500;      	  % pressure altitude in the stationary flight condition [m]
V0     = 97;          % true airspeed in the stationary flight condition [m/sec]
beta0 = 0.1;
phi0 = 0.1;
th0    = 0.0;       	  % pitch angle in the stationary flight condition [rad]
m      = 6000;           % mass [kg]
da_input = 0.1/180*pi;
dr_input = 0.1/180*pi;
d_time = 200;
t = 200;

% calc
asymmetric = @ss_a;
[sys_a3,eig_a3] = ss_a(V0,hp0,m,rho0,lambda,Temp0,g,R,S,b,CYbdot,KX2,KXZ,KZ2,Cnbdot,CYb,CYp,CYr,Clb,Clp,Clr,Cnb,Cnp,Cnr,CYda,CYdr,Clda,Cldr,Cnda,Cndr,C,D,th0);

dt5 = [1:d_time/dt];
t5 = 0:dt:t;
u5 = da_input(ones(length(t5),2));
u5(dt5,1) = da_input;
u5(dt5,2) = dr_input;
% u2 = de_input(ones(de_time/dt,1));
% u2 = [u2; zeros(t/dt-length(u2)+1,1)];
x0 = [beta0; phi0; 0; 0];
figure(5);
lsim(sys_a3,u5,t5,x0)

%% Functions
function [sys_s,eig_symmetric] = ss_s(V0,hp0,m,rho0,lambda,Temp0,g,R,S,c,CZadot,Cmadot,KY2,CXu,CXa,CXq,CZu,CZa,CZq,Cmu,Cma,Cmq,CXde,CZde,Cmde,C,D,th0)
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

function [sys_a,eig_asymmetric] = ss_a(V0,hp0,m,rho0,lambda,Temp0,g,R,S,b,CYbdot,KX2,KXZ,KZ2,Cnbdot,CYb,CYp,CYr,Clb,Clp,Clr,Cnb,Cnp,Cnr,CYda,CYdr,Clda,Cldr,Cnda,Cndr,C,D,th0)
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
    sys_a = ss(A1_asymmetric,B1_asymmetric,C,D);
end

% % Short period
% figure(1);
% subplot(2,2,1);
% plot(t1,y(:,1))                                        % V against time
% 
% subplot(2,2,2);
% plot(t1,y(:,2))                                        % AOA against time
% 
% subplot(2,2,3);
% plot(t1,y(:,3))                                        % Pitch angle against time
% 
% subplot(2,2,4);
% plot(t1,y(:,4))                                        % Pitch rate against time

% % symmetric_eq1 = [CXu-2*muc*c/V0*(ddt), CXa, CZ0, CXq; CZu, CZa+(CZadot-2*muc)*c/V0*(ddt), -CX0, CZq+2*muc; 0, 0, -c/V0*(ddt), 1; Cmu, Cma+Cmadot*c/V0*(ddt),0,Cmq-2*muc*KY2]
% % symmetric_eq2 = [ucirc;alpha;theta;q*c/V0]
% % sym1 = symmetric_eq1*symmetric_eq2

