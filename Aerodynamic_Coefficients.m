% -------------------------------- Run Matlab Files
T1 = excel_data_reader.T1;
T2 = excel_data_reader.T2;
T3 = excel_data_reader.T3;
Thrustdata = importdata('thrust.dat');

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
% % Basic Empty Weight [kg]
% m_BEM = 9165.0*0.453592; 
% % Payload Weight [kg]
% m_payload = 695; %REFERENCE DATA MASS
% % Block Fuel [kg]
% m_block_fuel = 4050*0.453592;
% % Mass Fuel Used [kg]
% m_fuel_used = T1(:,9)*0.453592;
% m_test = m_BEM+m_payload+m_block_fuel-m_fuel_used;
% Weight dependent on time
W = (massbalance.Weight1)*g; % CHANGE TO THE MASS GIVEN BY ROWAN
% Stick Force [N]
F_e = T2(:,9);
% Dimensionless Thrust Moment Arm 
CmTc = - 0.0064;
% Radius Engine [m]
Radius_Engine = 34.29/100;






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
CL_alpha = (Fit_CL_alpha(1))*(180/pi); % [1/rad]
alpha0 = -Fit_CL_alpha(2)/(Fit_CL_alpha(1)); % [deg]

%---------------------------------- CD-alpha curve

Thrust = sum(Thrustdata(1:6,1:2),2);
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
C_D_corr = CD0 + ((CL_alpha/(180/pi))*(alpha_array-alpha0)).^2/(pi*A*e);

%----------------------- Range
M_end = round(M(end),2,'significant');
M_first = round(M(1),2,'significant');
M_Range = [M_end M_first];

mu_ref = 1.716*10^-5 ;
S_suther = 110.4;
C1 = 1.458*10^-6;

mu = (C1*T.^(1.5))./(T+S_suther);
Re = (rho.*V_t*c)./mu;

Re_end = Re(end);
Re_first = Re(1);

Re_Range = [Re_end Re_first];

%------------------------ Stationary Measurements Elevator Trim Curve
%Everything that ends on _ET = Elevator Trim

% Pressure altitude in the stationary flight condition [m]
hp0_ET    = T2(:,4)*0.3048;
% Calibrated airspeed in stationary flight conidition [m/s]
Vc_ET     = (T2(:,5)-2)*0.514444; 
% Total Temperature [K]
T_m_ET    = T2(:,13)+273.15; 
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
%Reduced Equivalent Airspeed V_e(Tilde)
V_e_red = V_e_ET.*sqrt(W_s./W_ET);
%Calculate reduced elevator deflection

%---------------------------------- Shift in center of gravity (SCG), OUTPUTS: Cm_delta &
%Cm_alpha
% Pressure altitude in the stationary flight condition [m]
hp0_SCG    = T3(:,4)*0.3048;
% Calibrated airspeed in stationary flight conidition [m/s]
Vc_SCG     = (T3(:,5)-2)*0.514444; 
% Total Temperature [K]
T_m_SCG    = T3(:,13)+273.15;
% Angle of Attack [deg]
alpha_array_SCG = T3(:,6); 
% Weight dependent on time
W_SCG = (massbalance.Weight3)*g; % CHANGE TO THE MASS GIVEN BY ROWAN
% Measured elevator deflection
delta_e_SCG = T3(:,7);

% Pressure Calculation
p_SCG = p0*(1.0 + lambda*hp0_SCG/Temp0).^(-g/(lambda*R)); %pressure [Pa]
% Density Calculation
rho_SCG    = rho0*((1+(lambda*hp0_SCG/Temp0))).^(-((g/(lambda*R))+1));   % [kg/m^3]  (air density)
% Mach number
M_SCG = sqrt(2.0/(gamma-1.0)*((1.0+(p0./p_SCG).*((1.0 + (gamma-1.0)/(2.0*gamma)*(rho0/p0*Vc_SCG.^2.0)).^(gamma/(gamma-1.0))-1.0)).^((gamma-1.0)/gamma)-1.0));
%Static Temperature
T_SCG = T_m_SCG./(1.0+((gamma-1.0)/2.0)*M_SCG.^2.0);
%Speed of Sound
a_SCG = sqrt(gamma*R*T_SCG);
%True Airspeed
V_t_SCG = a_SCG.*M_SCG;
%Normal Coefficient C_N. delta C_N is neglected
C_N = W_SCG(1)/(0.5*rho_SCG(1)*V_t_SCG(1)^2*S);
C_m_delta = -(1.0 /((delta_e_SCG(2)-delta_e_SCG(1))*0.0174533))*C_N*(((massbalance.xcg3(1)-massbalance.xcg3(2)))/c);


%---------------------------------- Reduced Elevator Control Force Curve
%Reduced elevator control force
F_e_red = F_e.*(W_s./(W_ET));

%---------------------------------- Elevator Trime Curve CONTINUED
T_ET_inter = sum(Thrustdata,2);
T_ET = T_ET_inter(7:end-((size(Thrustdata,1)/2)+2));


Ts2_inter = Thrustdata(:,1);
Ts2 = Ts2_inter(((size(Thrustdata,1)/2)+7):end-2);


%Dimensionless thrust coefficient T_c
Tc2 = T_ET./(0.5*rho_ET.*V_t_ET.^2*pi*Radius_Engine^2);
%Dimensionless standard thrust coefficient T_cs
Tcs2 = 2*Ts2./(0.5*rho_ET.*V_t_ET.^2*pi*Radius_Engine^2);

delta_e_red = delta_e*0.0174533 - (1/C_m_delta)*CmTc*(Tcs2-Tc2);

%---------------------------------- Aanpassen Volgorde van lage snelheid
%naar hoge snelheid, DEZE AANPASSEN ALS JE DE FLIGHT DATA INPUT

% Eerst een table maken die bestaat uit de volgende parameters:
% | V_e_red [m/s] | delta_e_red [rad] | F_e_red [N] | alpha_array_ET [deg]

Variables_Aero_Coeff = [V_e_red delta_e_red F_e_red alpha_array_ET*0.0174533];

Final_Variables_Aero_Coeff = [
    Variables_Aero_Coeff(4,:),
    Variables_Aero_Coeff(3,:),
    Variables_Aero_Coeff(2,:),
    Variables_Aero_Coeff(1,:),
    Variables_Aero_Coeff(5,:),
    Variables_Aero_Coeff(6,:)
    ];
%,Variables_Aero_Coeff(7,:)

%---------------------------------- Longitudinal Stability, C_m_alpha
% Assumption: Lineaire functie 
Fit_C_m_alpha = polyfit(Final_Variables_Aero_Coeff(:,4),Final_Variables_Aero_Coeff(:,2),1);

C_m_alpha = Fit_C_m_alpha(1)*-C_m_delta;

%---------------------------------- Plot graphs
% % C_L - Alpha Curve
% figure
% hold on
% 
% subplot(2,2,1)
% plot(alpha_array,C_L,'b',x,yfit,'g')
% ax = gca;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
% legend('C_L-\alpha','C_L-\alpha extrapolated')
% xlabel('\alpha [deg]')
% ylabel('C_L [-]')
% title(['C_L-\alpha, M = ',num2str(M_Range(1)),' - ',num2str(M_Range(2)),' Re = ',num2str(Re_Range(1)),'-',])
% 
% % C_L - C_D Curve
% subplot(2,2,2)
% plot(C_D,C_L,'b',C_D_corr,C_L,'r')
% ax = gca;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
% legend('C_L-C_D','C_L-C_D_{corr}')
% xlabel('C_D [-]')
% ylabel('C_L [-]')
% title(['C_L-C_D, M = ',num2str(M_Range(1)),'-',num2str(M_Range(end)),' Re = ',num2str(Re_Range(1)),'-',num2str(Re_Range(end))]) 
% 
% % CL^2-CD Curve
% subplot(2,2,3)
% plot(C_L2,C_D,'b',x2,yfit2,'g')
% ax = gca;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
% legend('C_L^2-C_D','C_L^2-C_D extrapolated')
% xlabel('C_L^2 [-]')
% ylabel('C_D [-]')
% title('C_D-C_L^2 plot')
% 
% % C_D - Alpha Curve
% subplot(2,2,4)
% plot(alpha_array,C_D,'b',alpha_array,C_D_corr,'r')
% ax = gca;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
% legend('C_D-\alpha','C_D_{corr}-\alpha')
% xlabel('\alpha [deg]')
% ylabel('C_D [-]')
% title('C_D-\alpha plot')

% Reduced Elevator Control Force Curve, V_e_red vs. F_e_red
figure
plot(Final_Variables_Aero_Coeff(:,1),Final_Variables_Aero_Coeff(:,3))
ax = gca;
ax.YDir = 'reverse';
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
legend('F_{e_{red}}-V_{e_{red}}')
xlabel('V_{e_{red}} [m/s]')
ylabel('F_{e_{red}} (-) [N]')
title('Reduced Elevator Control Force')


% Reduced Elevator Trim Curve, V_e_red vs. delta_e_red
figure
plot(Final_Variables_Aero_Coeff(:,1),Final_Variables_Aero_Coeff(:,2))
ax = gca;
ax.YDir = 'reverse';
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
legend('\delta_{e_{red}}-V_{e_{red}}')
xlabel('V_{e_{red}} [m/s]')
ylabel('\delta_{e_{red}} (-) [rad]')
%title('Reduced Elevator Trime Curve')

figure
% Reduced Elevator, alpha vs. delta_e_red
plot(Final_Variables_Aero_Coeff(:,4),Final_Variables_Aero_Coeff(:,2))
ax = gca;
ax.YDir = 'reverse';
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
legend('\delta_{e_{red}}-\alpha')
xlabel('\alpha [rad]')
ylabel('\delta_{e_{red}} (-) [rad]')
%title('Reduced Elevator Trime Curve alpha vs delta e red')

%CL alpha
figure 
plot(alpha_array,C_L,'b',x,yfit,'r')
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
legend('C_L-\alpha','C_L-\alpha extrapolated')
xlabel('\alpha [deg]')
ylabel('C_L [-]')

%CL CD
figure
plot(C_D,C_L,'b',C_D_corr,C_L,'r')
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
legend('C_L-C_D','C_L-C_D_{corr}')
xlabel('C_D [-]')
ylabel('C_L [-]')

CLa = CL_alpha;
Cma = C_m_alpha;
Cmde = C_m_delta;

%clearvars -except alpha0 C_m_alpha C_m_delta CL_alpha e
clearvars a a_ET a_SCG alpha0 slpha_array alpha_array_ET ...
    alpha_SCG ands bh blockfuel C1 C_D C_D_corr C_L C_L2 ...
    C_m_alpha C_m_delta C_N CD0 CL-alpha CmTc delta_e delta_e_red ...
    delta_e_SCG e excel_data_reader F_e F_e_red Fit_C_m_alpha Fit_CL2_CD ...
    Fit_CL-alpha gamma hp0 hp0_ET hp0_SCG ih lh lh_c M M_end ...
    M_ET M_first M_Range M_SCG mu mu_ref p p0 p_ET p_SCG Radius_Engine ...
    Re Re_end Re_first Re_Range rho rho_ET rho_SCG S_suther Sh ...
    Sh_S T T1 T2 T3 T_ET T_ET_inter T_m T_m_ET T_m_SCG T_SCG Tc2 Tcs2 ...
    Thrust Thrustdata Ts2 Ts2_inter V_e V_e_ET V_e_red V_t ...
    V_t_ET V_t_SCG Vc Vc_ET Vc_SCG Vh_V W W_ET W_s W_SCG x x2 yfit yfit2