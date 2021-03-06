% Citation 550 - Linear simulation

% xcg = 0.25*c
%Obtain inputs
T1 = excel_data_reader.T1;
T2 = excel_data_reader.T2;
T3 = excel_data_reader.T3;

% Stationary flight condition
hp0    = 1; 
hp1    = T1(:,4)*0.3048;      	  % pressure altitude in the stationary flight condition [m]
hp2    = T2(:,4)*0.3048; 
hp3    = T3(:,4)*0.3048;
V0     = 1;            % true airspeed in the stationary flight condition [m/sec]
alpha0 = 1;       	  % angle of attack in the stationary flight condition [rad]
th0    = 1;       	  % pitch angle in the stationary flight condition [rad]

Vc1     = (T1(:,5)-2)*0.514444; %calibrated airspeed in stationary flight conidition [m/s]
Vc2     = (T2(:,5)-2)*0.514444;
Vc3     = (T3(:,5)-2)*0.514444;
% Aircraft mass
m      = 1;         	  % mass [kg]

% aerodynamic properties
e      = 1;            % Oswald factor [ ]
CD0    = 1;            % Zero lift drag coefficient [ ]
CLa    = 1;            % Slope of CL-alpha curve [ ]

% Longitudinal stability
Cma    = 1;            % longitudinal stabilty [ ]
Cmde   = 1;            % elevator effectiveness [ ]

% Aircraft geometry

S      = 30.00;	          % wing area [m^2]
Sh     = 0.2*S;           % stabiliser area [m^2]
Sh_S   = Sh/S;	          % [ ]
lh     = 0.71*5.968;      % tail length [m]
c      = 2.0569;	  % mean aerodynamic cord [m]
lh_c   = lh/c;	          % [ ]
b      = 15.911;	  % wing span [m]
bh     = 5.791;	          % stabilser span [m]
A      = b^2/S;           % wing aspect ratio [ ]
Ah     = bh^2/Sh;         % stabilser aspect ratio [ ]
Vh_V   = 1;		  % [ ]
ih     = -2*pi/180;       % stabiliser angle of incidence [rad]

% Constant values concerning atmosphere and gravity

rho0   = 1.2250;          % air density at sea level [kg/m^3] 
lambda = -0.0065;         % temperature gradient in ISA [K/m]
Temp0  = 288.15;          % temperature at sea level in ISA [K]
R      = 287.05;          % specific gas constant [m^2/sec^2K]
g      = 9.81;            % [m/sec^2] (gravity constant)
gamma  = 1.4;             % air ratio specific heats 
p0     = 101325;          % sea level pressure [Pa]
T_m1    = T1(:,10)+273.15;
T_m2    = T2(:,10)+273.15;
T_m3    = T3(:,10)+273.15;

rho    = rho0*((1+(lambda*hp0/Temp0)))^(-((g/(lambda*R))+1));   % [kg/m^3]  (air density)
W      = m*g;				                        % [N]       (aircraft weight)

Ws     = 60500; % [N]


% Pressure calculations
p1 = p0*(1.0 + lambda*hp1/Temp0).^(-g/(lambda*R)); %pressure [Pa]
p2 = p0*(1.0 + lambda*hp2/Temp0).^(-g/(lambda*R));
p3 = p0*(1.0 + lambda*hp3/Temp0).^(-g/(lambda*R));

% Density calulations 
rho1    = rho0*((1+(lambda*hp1/Temp0))).^(-((g/(lambda*R))+1));   % [kg/m^3]  (air density)
rho2    = rho0*((1+(lambda*hp2/Temp0))).^(-((g/(lambda*R))+1));
rho3    = rho0*((1+(lambda*hp3/Temp0))).^(-((g/(lambda*R))+1));

% Mach number
A1 = (p0./p1);
Mach.M1 = sqrt(2.0/(gamma-1.0)*((1.0+A1.*((1.0 + (gamma-1.0)/(2.0*gamma)*(rho0/p0*Vc1.^2.0)).^(gamma/(gamma-1.0))-1.0)).^((gamma-1.0)/gamma)-1.0));

A2 = (p0./p2);
Mach.M2 = sqrt(2.0/(gamma-1.0)*((1.0+A2.*((1.0 + (gamma-1.0)/(2.0*gamma)*(rho0/p0*Vc2.^2.0)).^(gamma/(gamma-1.0))-1.0)).^((gamma-1.0)/gamma)-1.0));

A3 = (p0./p3);
Mach.M3 = sqrt(2.0/(gamma-1.0)*((1.0+A3.*((1.0 + (gamma-1.0)/(2.0*gamma)*(rho0/p0*Vc3.^2.0)).^(gamma/(gamma-1.0))-1.0)).^((gamma-1.0)/gamma)-1.0));

clearvars T1 T2 T3 hp0 hp1 hp2 hp3 V0 alpha0 th0 Vc1 Vc2 Vc3 m e CD0 CLa ...
          Cma Cmde S Sh Sh_S lh c lh_c b bh A Ah Vh_V ih rho0 lambda ...
          Temp0 R g gamma p0 T_m1 T_m2 T_m3 rho W Ws p1 p2 p3 rho1 ...
          rho2 rho3 A1 A2 A3 a1 a2 a3 V_t1 V_t2 V_t3 V_e1 V_e2 V_e3 ...
          Tabel1 Tabel2 Tabel3
