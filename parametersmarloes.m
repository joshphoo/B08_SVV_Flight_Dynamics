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
parameter.M1 = sqrt(2.0/(gamma-1.0)*((1.0+A1.*((1.0 + (gamma-1.0)/(2.0*gamma)*(rho0/p0*Vc1.^2.0)).^(gamma/(gamma-1.0))-1.0)).^((gamma-1.0)/gamma)-1.0));

A2 = (p0./p2);
parameter.M2 = sqrt(2.0/(gamma-1.0)*((1.0+A2.*((1.0 + (gamma-1.0)/(2.0*gamma)*(rho0/p0*Vc2.^2.0)).^(gamma/(gamma-1.0))-1.0)).^((gamma-1.0)/gamma)-1.0));

A3 = (p0./p3);
parameter.M3 = sqrt(2.0/(gamma-1.0)*((1.0+A3.*((1.0 + (gamma-1.0)/(2.0*gamma)*(rho0/p0*Vc3.^2.0)).^(gamma/(gamma-1.0))-1.0)).^((gamma-1.0)/gamma)-1.0));

%Static Temperature
T1 = T_m1./(1.0+((gamma-1.0)/2.0)*parameter.M1.^2.0);
T2 = T_m2./(1.0+((gamma-1.0)/2.0)*parameter.M2.^2.0);
T3 = T_m3./(1.0+((gamma-1.0)/2.0)*parameter.M3.^2.0);

%Speed of Sound
a1 = sqrt(gamma*R*T1);
a2 = sqrt(gamma*R*T2);
a3 = sqrt(gamma*R*T3);

%True Airspeed
V_t1 = a1.*parameter.M1;
V_t2 = a2.*parameter.M2;
V_t3 = a3.*parameter.M3;

%Equivalent Airspeed
V_e1 = V_t1.*sqrt(rho1./rho0);
V_e2 = V_t2.*sqrt(rho2./rho0);
V_e3 = V_t3.*sqrt(rho3./rho0);

run('MassBalance.m'); 
%Reduced equivalent airspeed  
V_re1 = V_e1.*sqrt(massbalance.Weight1./Ws); 
V_re2 = V_e2.*sqrt(massbalance.Weight2./Ws); 
V_re3 = V_e3.*sqrt(massbalance.Weight3./Ws); 

R = 34.29/100;

Thrustdata = importdata('thrust.dat');

Tt1 = [sum(Thrustdata(1,:)),
    sum(Thrustdata(2,:)),
    sum(Thrustdata(3,:)),
    sum(Thrustdata(4,:)),
    sum(Thrustdata(5,:)),
    sum(Thrustdata(6,:))];

Tt2 = [sum(Thrustdata(7,:)),
    sum(Thrustdata(8,:)),
    sum(Thrustdata(9,:)),
    sum(Thrustdata(10,:)),
    sum(Thrustdata(11,:)),
    sum(Thrustdata(12,:)),
    sum(Thrustdata(13,:))];

Tt3 = [sum(Thrustdata(14,:)),
    sum(Thrustdata(15,:))];

Ts1 = [Thrustdata(16),
    Thrustdata(17),
    Thrustdata(18),
    Thrustdata(19),
    Thrustdata(20),
    Thrustdata(21)];

Ts2 = [Thrustdata(22),
    Thrustdata(23),
    Thrustdata(24),
    Thrustdata(25),
    Thrustdata(26),
    Thrustdata(27),
    Thrustdata(28)];

Ts3 = [Thrustdata(29),
    Thrustdata(30)];

%Dimensionless thrust coefficient T_c
Tc1 = Tt1./(0.5*rho1.*V_t1.^2*pi*R^2);
Tc2 = Tt2./(0.5*rho2.*V_t2.^2*pi*R^2);
Tc3 = Tt3./(0.5*rho3.*V_t3.^2*pi*R^2);

%Dimensionless standard thrust coefficient T_cs
Tcs1 = 2*Ts1./(0.5*rho1.*V_t1.^2*pi*R^2);
Tcs2 = 2*Ts2./(0.5*rho2.*V_t2.^2*pi*R^2);
Tcs3 = 2*Ts3./(0.5*rho3.*V_t3.^2*pi*R^2);

Cmdelta = -1.1642;
CmTc = - 0.0064
Tabel1 = excel_data_reader.T1;
Tabel2 = excel_data_reader.T2;
Tabel3 = excel_data_reader.T3;
deltaeeqmeas2 = Tabel2(:,7);
deltaeeqmeas3 = Tabel3(:,7);
%reduced elevator deflection 
%deltaeeq1 = deltaeeqmeas - 1/Cmdelta*CmTc*(Tcs1-Tc1);
deltaeeq2 = 0.0174533*deltaeeqmeas2 - 1/Cmdelta*CmTc*(Tcs2-Tc2);
deltaeeq3 = 0.0174533*deltaeeqmeas3 - 1/Cmdelta*CmTc*(Tcs3-Tc3);

plot(deltaeeq2,V_re2)
clearvars T1 T2 T3 hp0 hp1 hp2 hp3 V0 alpha0 th0 Vc1 Vc2 Vc3 m e CD0 CLa ...
          Cma Cmde S Sh Sh_S lh c lh_c b bh A Ah Vh_V ih rho0 lambda ...
          Temp0 R g gamma p0 T_m1 T_m2 T_m3 rho W Ws p1 p2 p3 rho1 ...
          rho2 rho3 A1 A2 A3 a1 a2 a3 V_t1 V_t2 V_t3 V_e1 V_e2 V_e3 ...
          Tabel1 Tabel2 Tabel3
