% close all;
% clear;
load('matlab.mat', 'flightdata')

motion.t1 = 0*3600 + 53*60 + 57;
motion.t2 = 1*3600 + 0*60 + 35;
motion.t3 = 1*3600 + 1*60 + 57;
motion.t4 = 1*3600 + 2*60 + 47;
motion.t5 = 0*3600 + 59*60 + 10;
motion.t6 = 1*3600 + 5*60 + 20;
motion.idx1 = find(flightdata.time.data==motion.t1);               % Starting time
motion.idxe1 = find(flightdata.time.data==motion.t1+220); 
motion.idx2 = find(flightdata.time.data==motion.t2);            % Starting time
motion.idxe2 = find(flightdata.time.data==motion.t2+50);  
motion.idx3 = find(flightdata.time.data==motion.t3);               % Starting time
motion.idxe3 = find(flightdata.time.data==motion.t3+45);
motion.idx4 = find(flightdata.time.data==motion.t4);               % Starting time
motion.idxe4 = find(flightdata.time.data==motion.t4+45);
motion.idx5 = find(flightdata.time.data==motion.t5);               % Starting time
motion.idxe5 = find(flightdata.time.data==motion.t5+50);
motion.idx6 = find(flightdata.time.data==motion.t6);               % Starting time
motion.idxe6 = find(flightdata.time.data==motion.t6+220);    
m_max = (9165+4050)*0.453592+95+92+74+66+61+75+78+86+68;


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
hp0    = flightdata.Dadc1_alt.data(motion.idx1)*0.3048;      	  % pressure altitude in the stationary flight condition [m]
V0     = flightdata.Dadc1_tas.data(motion.idx1)*0.514444;          % true airspeed in the stationary flight condition [m/sec]
%flightdata.vane_AOA.data(motion.idx1)/180*pi();    % angle of attack in the stationary flight condition [rad]
th0    = flightdata.Ahrs1_Pitch.data(motion.idx1)/180*pi();       	  % pitch angle in the stationary flight condition [rad]
m      = m_max-(flightdata.rh_engine_FU.data(motion.idx1)+flightdata.lh_engine_FU.data(motion.idx1))*0.453592;           % mass [kg]

t = flightdata.time.data(motion.idx1:motion.idxe1);
t = t-t(1);
de = flightdata.delta_e.data(motion.idx1:motion.idxe1)/180*pi();
x0 = [1; alpha0; th0; 0];

symmetric = @ss_s;
[sys_s1,eig_s1] = symmetric(V0,hp0,m,rho0,lambda,Temp0,g,R,S,c,CZadot,Cmadot,KY2,CXu,CXa,CXq,CZu,CZa,CZq,Cmu,Cma,Cmq,CXde,CZde,Cmde,C,D,th0);

figure(1);
lsim(sys_s1,de,t,x0)

%% Phugoid
% input
hp0    = flightdata.Dadc1_alt.data(motion.idx2);      	  % pressure altitude in the stationary flight condition [m]
V0     = flightdata.Dadc1_tas.data(motion.idx2);          % true airspeed in the stationary flight condition [m/sec]
alpha0 = flightdata.vane_AOA.data(motion.idx2);    % angle of attack in the stationary flight condition [rad]
th0    = flightdata.Ahrs1_Pitch.data(motion.idx2)/180*pi();       	  % pitch angle in the stationary flight condition [rad]
m      = m_max-(flightdata.rh_engine_FU.data(motion.idx2)+flightdata.lh_engine_FU.data(motion.idx2))*0.453592;           % mass [kg]

t = flightdata.time.data(motion.idx2:motion.idxe2);
t = t-t(1);
de = flightdata.delta_e.data(motion.idx2:motion.idxe2)/180*pi();
x0 = [1; alpha0; th0; 0];

symmetric = @ss_s;
[sys_s2,eig_s2] = symmetric(V0,hp0,m,rho0,lambda,Temp0,g,R,S,c,CZadot,Cmadot,KY2,CXu,CXa,CXq,CZu,CZa,CZq,Cmu,Cma,Cmq,CXde,CZde,Cmde,C,D,th0);

figure(2);
lsim(sys_s2,de,t,x0)


%% Dutch Roll
% input
hp0    = flightdata.Dadc1_alt.data(motion.idx3);      	  % pressure altitude in the stationary flight condition [m]
V0     = flightdata.Dadc1_tas.data(motion.idx3);          % true airspeed in the stationary flight condition [m/sec]
m      = m_max-(flightdata.rh_engine_FU.data(motion.idx3)+flightdata.lh_engine_FU.data(motion.idx3))*0.453592;           % mass [kg]
beta0 = 0;
phi0 = flightdata.Ahrs1_Roll.data(motion.idx3)/180*pi();

t = flightdata.time.data(motion.idx3:motion.idxe3);
t = t-t(1);
dar = [];
dar(:,1) = flightdata.delta_a.data(motion.idx3:motion.idxe3)/180*pi();
dar(:,2) = flightdata.delta_r.data(motion.idx3:motion.idxe3)/180*pi();
x0 = [beta0; phi0; 0; 0];

% calc
asymmetric = @ss_a;
[sys_a1,eig_a1] = ss_a(V0,hp0,m,rho0,lambda,Temp0,g,R,S,b,CYbdot,KX2,KXZ,KZ2,Cnbdot,CYb,CYp,CYr,Clb,Clp,Clr,Cnb,Cnp,Cnr,CYda,CYdr,Clda,Cldr,Cnda,Cndr,C,D);

figure(3);
lsim(sys_a1,dar,t,x0)

%% Dutch Roll Damped
% input
hp0    = flightdata.Dadc1_alt.data(motion.idx4);      	  % pressure altitude in the stationary flight condition [m]
V0     = flightdata.Dadc1_tas.data(motion.idx4);          % true airspeed in the stationary flight condition [m/sec]
m      = m_max-(flightdata.rh_engine_FU.data(motion.idx4)+flightdata.lh_engine_FU.data(motion.idx4))*0.453592;           % mass [kg]
beta0 = 0;
phi0 = flightdata.Ahrs1_Roll.data(motion.idx4)/180*pi();

t = flightdata.time.data(motion.idx4:motion.idxe4);
t = t-t(1);
dar = [];
dar(:,1) = flightdata.delta_a.data(motion.idx4:motion.idxe4)/180*pi();
dar(:,2) = flightdata.delta_r.data(motion.idx4:motion.idxe4)/180*pi();
x0 = [beta0; phi0; 0; 0];

% calc
asymmetric = @ss_a;
[sys_a2,eig_a2] = ss_a(V0,hp0,m,rho0,lambda,Temp0,g,R,S,b,CYbdot,KX2,KXZ,KZ2,Cnbdot,CYb,CYp,CYr,Clb,Clp,Clr,Cnb,Cnp,Cnr,CYda,CYdr,Clda,Cldr,Cnda,Cndr,C,D);

figure(4);
lsim(sys_a2,dar,t,x0)

%% Aperiodic Roll
% input
hp0    = flightdata.Dadc1_alt.data(motion.idx5);      	  % pressure altitude in the stationary flight condition [m]
V0     = flightdata.Dadc1_tas.data(motion.idx5);          % true airspeed in the stationary flight condition [m/sec]
m      = m_max-(flightdata.rh_engine_FU.data(motion.idx5)+flightdata.lh_engine_FU.data(motion.idx5))*0.453592;           % mass [kg]
beta0 = 0;
phi0 = flightdata.Ahrs1_Roll.data(motion.idx5)/180*pi();          % mass [kg]

t = flightdata.time.data(motion.idx5:motion.idxe5);
t = t-t(1);
dar = [];
dar(:,1) = flightdata.delta_a.data(motion.idx5:motion.idxe5)/180*pi();
dar(:,2) = flightdata.delta_r.data(motion.idx5:motion.idxe5)/180*pi();
x0 = [beta0; phi0; 0; 0];

% calc
asymmetric = @ss_a;
[sys_a3,eig_a3] = ss_a(V0,hp0,m,rho0,lambda,Temp0,g,R,S,b,CYbdot,KX2,KXZ,KZ2,Cnbdot,CYb,CYp,CYr,Clb,Clp,Clr,Cnb,Cnp,Cnr,CYda,CYdr,Clda,Cldr,Cnda,Cndr,C,D);

figure(5);
lsim(sys_a3,dar,t,x0)

%% Spiral
% input
hp0    = flightdata.Dadc1_alt.data(motion.idx6);      	  % pressure altitude in the stationary flight condition [m]
V0     = flightdata.Dadc1_tas.data(motion.idx6);          % true airspeed in the stationary flight condition [m/sec]
m      = m_max-(flightdata.rh_engine_FU.data(motion.idx6)+flightdata.lh_engine_FU.data(motion.idx6))*0.453592;           % mass [kg]
beta0 = 0;
phi0 = flightdata.Ahrs1_Roll.data(motion.idx6)/180*pi(); 

t = flightdata.time.data(motion.idx6:motion.idxe6);
t = t-t(1);
dar = [];
dar(:,1) = flightdata.delta_a.data(motion.idx6:motion.idxe6)/180*pi();
dar(:,2) = flightdata.delta_r.data(motion.idx6:motion.idxe6)/180*pi();
x0 = [beta0; phi0; 0; 0];

% calc
asymmetric = @ss_a;
[sys_a4,eig_a4] = ss_a(V0,hp0,m,rho0,lambda,Temp0,g,R,S,b,CYbdot,KX2,KXZ,KZ2,Cnbdot,CYb,CYp,CYr,Clb,Clp,Clr,Cnb,Cnp,Cnr,CYda,CYdr,Clda,Cldr,Cnda,Cndr,C,D);

figure(6);
lsim(sys_a4,dar,t,x0)


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

