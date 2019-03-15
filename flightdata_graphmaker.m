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
%Motion 1
idx1 = find(flightdata.time.data==9);               % Starting time
figure(1);
subplot(2,2,1);
xs = flightdata.time.data(idx1:end);
ys1 = flightdata.Dadc1_mach.data(idx1:end);
plot(xs,ys1)                                        % Mach against time

subplot(2,2,2); 
ys2 = flightdata.vane_AOA.data(idx1:end);
plot(xs,ys2)                                        % AOA against time

subplot(2,2,3)
ys3 = flightdata.Ahrs1_Pitch.data(idx1:end);
plot(xs,ys3)                                        % Pitch against time

subplot(2,2,4)
ys4 = flightdata.Ahrs1_bPitchRate.data(idx1:end);
plot(xs,ys4)                                        % Pitch rate against time

hold on

%Motion 2
idx2 = find(flightdata.time.data==1500);            % Starting time
figure(2);
subplot(2,2,1);
xs = flightdata.time.data(idx2:end);
ys1 = flightdata.Dadc1_mach.data(idx2:end);
plot(xs,ys1)                                        % Mach against time

subplot(2,2,2); 
ys2 = flightdata.vane_AOA.data(idx2:end);
plot(xs,ys2)                                        % AOA against time

subplot(2,2,3)
ys3 = flightdata.Ahrs1_Pitch.data(idx2:end);
plot(xs,ys3)                                        % Pitch against time

subplot(2,2,4)
ys4 = flightdata.Ahrs1_bPitchRate.data(idx2:end);
plot(xs,ys4)                                        % Pitch rate against time

hold on

%Assymetric
%Motion 3
idx3 = find(flightdata.time.data==9);               % Starting time
figure(3)
subplot(2,2,1)
xa = flightdata.time.data(idx3:end);
ya1 = flightdata.Ahrs1_Roll.data(idx3:end);
plot(xa,ya1)                                        % Roll against time

subplot(2,2,2)
ya2 = flightdata.Ahrs1_bRollRate.data(idx3:end);
plot(xa,ya2)                                        % Roll rate against time

subplot(2,2,3)
ya3 = flightdata.Fms1_trueHeading.data(idx3:end);             %incorrect, how to obtain yaw?
plot(xa,ya3)                                        % Yaw against time

subplot(2,2,4)
ya4 = flightdata.Ahrs1_bYawRate.data(idx3:end);
plot(xa,ya4)

hold off

%Motion 4
idx4 = find(flightdata.time.data==3900);               % Starting time
figure(4)
subplot(2,2,1)
xa = flightdata.time.data(idx4:end);
ya1 = flightdata.Ahrs1_Roll.data(idx4:end);
plot(xa,ya1)                                        % Roll against time

subplot(2,2,2)
ya2 = flightdata.Ahrs1_bRollRate.data(idx4:end);
plot(xa,ya2)                                        % Roll rate against time

subplot(2,2,3)
ya3 = flightdata.Fms1_trueHeading.data(idx4:end);             %incorrect, how to obtain yaw?
plot(xa,ya3)                                        % Yaw against time

subplot(2,2,4)
ya4 = flightdata.Ahrs1_bYawRate.data(idx4:end);
plot(xa,ya4)

hold off

%Motion 5
idx5 = find(flightdata.time.data==4000);               % Starting time
figure(5)
subplot(2,2,1)
xa = flightdata.time.data(idx5:end);
ya1 = flightdata.Ahrs1_Roll.data(idx5:end);
plot(xa,ya1)                                        % Roll against time

subplot(2,2,2)
ya2 = flightdata.Ahrs1_bRollRate.data(idx5:end);
plot(xa,ya2)                                        % Roll rate against time

subplot(2,2,3)
ya3 = flightdata.Fms1_trueHeading.data(idx5:end);             %incorrect, how to obtain yaw?
plot(xa,ya3)                                        % Yaw against time

subplot(2,2,4)
ya4 = flightdata.Ahrs1_bYawRate.data(idx5:end);
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
















