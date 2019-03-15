% -------------------------------- Run Matlab Files
run('excell_data_reader.m');
T1 = excel_data_reader.T1;
T2 = excel_data_reader.T2;
Thrustdata = importdata('thrust.dat');
run('MassBalance.m');

% -------------------------------- Parameters
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
W_s    = 60500; % [N]

% Constant values concerning atmosphere and gravity

rho0   = 1.2250;          % air density at sea level [kg/m^3] 
lambda = -0.0065;         % temperature gradient in ISA [K/m]
Temp0  = 288.15;          % temperature at sea level in ISA [K]
R      = 287.05;          % specific gas constant [m^2/sec^2K]
g      = 9.81;            % [m/sec^2] (gravity constant)
gamma  = 1.4;             % air ratio specific heats 
p0     = 101325;          % sea level pressure [Pa]

%--------------------------------- Calculate Parameters
% Pressure altitude in the stationary flight condition [m]
hp0    = T1(:,4)*0.3048; 
% Calibrated airspeed in stationary flight conidition [m/s]
Vc     = (T1(:,5)-2)*0.514444; 
% Total Temperature [K]
T_m    = T1(:,10)+273.15; 
% Angle of Attack [deg]
alpha_array = T1(:,6); 
% Basic Empty Weight [kg]
m_BEM = 9165.0*0.453592; 
% Payload Weight [kg]
m_payload = 695; %REFERENCE DATA MASS
% Block Fuel [kg]
m_block_fuel = 4050*0.453592;
% Mass Fuel Used [kg]
m_fuel_used = T1(:,9)*0.453592;
m_test = m_BEM+m_payload+m_block_fuel-m_fuel_used;
% Weight dependent on time
W = (massbalance.Weight1)*g; % CHANGE TO THE MASS GIVEN BY ROWAN










% ----------------------------------- First Stationary Measurements CL-CD
% Pressure Calculation
p = p0*(1.0 + lambda*hp0/Temp0).^(-g/(lambda*R)); %pressure [Pa]

% Density Calculation
rho    = rho0*((1+(lambda*hp0/Temp0))).^(-((g/(lambda*R))+1));   % [kg/m^3]  (air density)

% Mach number
M = sqrt(2.0/(gamma-1.0)*((1.0+(p0./p).*((1.0 + (gamma-1.0)/(2.0*gamma)*(rho0/p0*Vc.^2.0)).^(gamma/(gamma-1.0))-1.0)).^((gamma-1.0)/gamma)-1.0));

%Static Temperature
T = T_m./(1.0+((gamma-1.0)/2.0)*M.^2.0);

%Speed of Sound
a = sqrt(gamma*R*T);

%True Airspeed
V_t = a.*M;

%Equivalent Airspeed
V_e = V_t.*sqrt(rho./rho0);

%------------------------------- Calculate C_L vs alpha curve
C_L = W./(0.5*rho.*V_t.^2*S);
alpha_array = alpha_array.';
C_L = C_L.'; 
Fit_CL_alpha = polyfit(alpha_array,C_L,1);
x = -5:0.01:10;
yfit = Fit_CL_alpha(1)*x+Fit_CL_alpha(2);
CL_alpha = Fit_CL_alpha(1); % [1/deg]
alpha0 = -Fit_CL_alpha(2)/CL_alpha; % [deg]

%---------------------------------- CD-alpha curve
Thrust = [
    sum(Thrustdata(1,:)),
    sum(Thrustdata(2,:)),
    sum(Thrustdata(3,:)),
    sum(Thrustdata(4,:)),
    sum(Thrustdata(5,:)),
    sum(Thrustdata(6,:))
    ]; % N
C_D = Thrust./(0.5*rho.*V_t.^2*S);
C_D = C_D.';
C_L2 = C_L.^2; %CL^2
Fit_CL2_CD = polyfit(C_L2,C_D,1);
x2 = -0.2:0.01:1;
yfit2 = Fit_CL2_CD(1)*x2+Fit_CL2_CD(2);


%Oswaldo factor
e = 1/(pi*A*Fit_CL2_CD(1));
% Zero lift drag coefficient
CD0 = Fit_CL2_CD(2);

% Corrected C_D with input: CD0 & e
C_D_corr = CD0 + (CL_alpha*(alpha_array-alpha0)).^2/(pi*A*e);

%----------------------- Range
M_Range = [M(end) M(1)]






%------------------------ Stationary Measurements Elevator Trim Curve
%Everything that ends on _ET = Elevator Trim

% Pressure altitude in the stationary flight condition [m]
hp0_ET    = T2(:,4)*0.3048;
% Calibrated airspeed in stationary flight conidition [m/s]
Vc_ET     = (T2(:,5)-2)*0.514444; 
% Total Temperature [K]
T_m_ET    = T2(:,10)+273.15; 
% Angle of Attack [deg]
alpha_array_ET = T2(:,6); 
% Weight dependent on time
W_ET = (massbalance.Weight2)*g; % CHANGE TO THE MASS GIVEN BY ROWAN
% Measured elevator deflection
delta_e = T2(:,7);

% Pressure Calculation
p_ET = p0*(1.0 + lambda*hp0_ET/Temp0).^(-g/(lambda*R)); %pressure [Pa]
% Density Calculation
rho_ET    = rho0*((1+(lambda*hp0_ET/Temp0))).^(-((g/(lambda*R))+1));   % [kg/m^3]  (air density)
% Mach number
M_ET = sqrt(2.0/(gamma-1.0)*((1.0+(p0./p_ET).*((1.0 + (gamma-1.0)/(2.0*gamma)*(rho0/p0*Vc_ET.^2.0)).^(gamma/(gamma-1.0))-1.0)).^((gamma-1.0)/gamma)-1.0));
%Static Temperature
T_ET = T_m_ET./(1.0+((gamma-1.0)/2.0)*M_ET.^2.0);
%Speed of Sound
a_ET = sqrt(gamma*R*T_ET);
%True Airspeed
V_t_ET = a_ET.*M_ET;
%Equivalent Airspeed
V_e_ET = V_t_ET.*sqrt(rho_ET./rho0);
%Reduced Equivalent Airspeed
V_e_red = V_e_ET.*sqrt(W_s./W_ET);
%Calculate reduced elevator deflection






%---------------------------------- Plot graphs
% C_L - Alpha Curve
subplot(2,2,1)
plot(alpha_array,C_L,'b',x,yfit,'g')
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
legend('C_L-\alpha','C_L-\alpha extrapolated')
xlabel('\alpha [deg]')
ylabel('C_L [-]')
t = 2
title(['C_L-\alpha plot, Mach range: ',num2str(M_Range(1)),'-',num2str(M_Range(2)),' Reynolds nr. range:'])

% C_L - C_D Curve
subplot(2,2,2)
plot(C_D,C_L,'b',C_D_corr,C_L,'r')
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
legend('C_L-C_D','C_L-C_D_{corr}')
xlabel('C_D [-]')
ylabel('C_L [-]')
title('C_L-C_D plot') 

% CL^2-CD Curve
subplot(2,2,3)
plot(C_L2,C_D,'b',x2,yfit2,'g')
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
legend('C_L^2-C_D','C_L^2-C_D extrapolated')
xlabel('C_L^2 [-]')
ylabel('C_D [-]')
title('C_D-C_L^2 plot')

% C_D - Alpha Curve
subplot(2,2,4)
plot(alpha_array,C_D,'b',alpha_array,C_D_corr,'r')
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
legend('C_D-\alpha','C_D_{corr}-\alpha')
xlabel('\alpha [deg]')
ylabel('C_D [-]')
title('C_D-\alpha plot')
