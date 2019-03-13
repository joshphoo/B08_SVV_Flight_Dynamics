%-----------------------FLIGHT DATA READER---------------------------------
%-----------------------Load the flightdata--------------------------------
load('matlab.mat', 'flightdata')
% T = table(flightdata.Dadc1_alt.data, flightdata.Dadc1_mach.data, flightdata.Dadc1_tat.data - flightdata.Dadc1_sat.data, flightdata.lh_engine_FMF.data, flightdata.rh_engine_FMF.data) ;
%-----------------Assign variables to necessary data-----------------------
%Velocity
%V = flightdata.
%General
% t = flightdata.time.data;                           % time(s)
% 
% M = flightdata.Dadc1_mach.data;                     % Mach number(-)

%Symmetric
figure(1);
subplot(2,2,1);
xs = flightdata.time.data;
ys1 = flightdata.Dadc1_mach.data;
plot(xs,ys1)                                        % Mach against time

subplot(2,2,2); 
ys2 = flightdata.vane_AOA.data;
plot(xs,ys2)                                        % AOA against time

subplot(2,2,3)
ys3 = flightdata.Ahrs1_Pitch.data;
plot(xs,ys3)                                        % Pitch against time

subplot(2,2,4)
ys4 = flightdata.Ahrs1_bPitchRate.data;
plot(xs,ys4)                                        % Pitch rate against time

hold on
%Assymetric
figure(2)
subplot(2,2,1)
xa = flightdata.time.data;
ya1 = flightdata.Ahrs1_Roll.data;
plot(xa,ya1)                                        % Roll against time

subplot(2,2,2)
ya2 = flightdata.Ahrs1_bRollRate.data;
plot(xa,ya2)                                        % Roll rate against time

subplot(2,2,3)
ya3 = flightdata.Fms1_trueHeading.data;             %incorrect, how to obtain yaw?
plot(xa,ya3)                                        % Yaw against time

subplot(2,2,4)
ya4 = flightdata.Ahrs1_bYawRate.data;
plot(xa,ya4)

hold off

% Phi = flightdata.Ahrs1_Roll.data;                   % Roll angle(deg)
% 
% p = flightdata.Ahrs1_bRollRate.data;                % Roll rate (deg/s)
% 
% %Psi = flightdata. It seems like this has to be calculated from flight heading    % Yaw(deg)
% 
% r = flightdata.Ahrs1_bYawRate.data;                 % Yaw rate(deg/s)


clearvars xa xs ya1 ya2 ya3 ya4 ys1 ys2 ys3 ys4
















