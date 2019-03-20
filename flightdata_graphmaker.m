%-----------------------FLIGHT DATA READER---------------------------------
%-----------------------Load the flightdata--------------------------------
clf
close all;
load('matlab.mat', 'flightdata')
% T = table(flightdata.Dadc1_alt.data, flightdata.Dadc1_mach.data, flightdata.Dadc1_tat.data - flightdata.Dadc1_sat.data, flightdata.lh_engine_FMF.data, flightdata.rh_engine_FMF.data) ;
%-----------------Assign variables to necessary data-----------------------
%Velocity
%V = flightdata.
%General
% t = flightdata.time.data;                           % time(s)
% 
% M = flightdata.Dadc1_mach.data;                     % Mach number(-)
motion.t1 = 0*3600 + 53*60 + 57;
motion.t2 = 1*3600 + 0*60 + 35;
motion.t3 = 1*3600 + 1*60 + 57;
motion.t4 = 1*3600 + 2*60 + 47;
motion.t5 = 0*3600 + 59*60 + 10;
motion.t6 = 1*3600 + 5*60 + 20;
%Symmetric
%Motion 1
motion.idx1 = find(flightdata.time.data==motion.t1-20);               % Starting time
motion.idxe1 = find(flightdata.time.data==motion.t1+220);              %temporary end time
refmotion(1) = figure(1);
title('Phugoid')
subplot(2,2,1);
xs = flightdata.time.data(motion.idx1:motion.idxe1)-flightdata.time.data(motion.idx1);;
ys1 = 0.514444.*flightdata.Dadc1_tas.data(motion.idx1:motion.idxe1);
plot(xs,ys1)                                        % TAS against time
title('TAS(m/s)')

subplot(2,2,2); 
ys2 = (pi/180).*flightdata.vane_AOA.data(motion.idx1:motion.idxe1);
plot(xs,ys2)                                        % AOA against time
title('AOA')

subplot(2,2,3)
ys3 = (pi/180).*flightdata.Ahrs1_Pitch.data(motion.idx1:motion.idxe1);
plot(xs,ys3)                                        % Pitch against time
title('Pitch')

subplot(2,2,4)
ys4 = (pi/180).*flightdata.Ahrs1_bPitchRate.data(motion.idx1:motion.idxe1);
plot(xs,ys4)                                        % Pitch rate against time
title('Pitch rate')

hold on

%Motion 2
motion.idx2 = find(flightdata.time.data==motion.t2);            % Starting time
motion.idxe2 = find(flightdata.time.data==motion.t2+50);              %temporary end time
refmotion(2) = figure(2);
title('Short period')
subplot(2,2,1);
xs = flightdata.time.data(motion.idx2:motion.idxe2)-flightdata.time.data(motion.idx2);
ys1 = 0.514444.*flightdata.Dadc1_tas.data(motion.idx2:motion.idxe2);
plot(xs,ys1)                                        % TAS against time
title('TAS(m/s)')

subplot(2,2,2); 
ys2 = (pi/180).*flightdata.vane_AOA.data(motion.idx2:motion.idxe2);
plot(xs,ys2)                                        % AOA against time
title('AOA')

subplot(2,2,3)
ys3 = (pi/180).*flightdata.Ahrs1_Pitch.data(motion.idx2:motion.idxe2);
plot(xs,ys3)                                        % Pitch against time
title('Pitch')

subplot(2,2,4)
ys4 = (pi/180).*flightdata.Ahrs1_bPitchRate.data(motion.idx2:motion.idxe2);
plot(xs,ys4)                                        % Pitch rate against time
title('Pitch rate')


hold on

%Assymetric
%Motion 3
motion.idx3 = find(flightdata.time.data==motion.t3);               % Starting time
motion.idxe3 = find(flightdata.time.data==motion.t3+45);              %temporary end time
refmotion(3) = figure(3);
title('Dutch roll')
subplot(2,2,1)
xa = flightdata.time.data(motion.idx3:motion.idxe3)-flightdata.time.data(motion.idx3);
ya1 = (pi/180).*flightdata.Ahrs1_Roll.data(motion.idx3:motion.idxe3);
plot(xa,ya1)                                        % Roll against time
title('Roll')

subplot(2,2,2)
ya2 = (pi/180).*flightdata.Ahrs1_bRollRate.data(motion.idx3:motion.idxe3);
plot(xa,ya2)                                        % Roll rate against time
title('Roll rate')

subplot(2,2,3)
ya3 = (pi/180).*flightdata.Ahrs1_bYawRate.data(motion.idx3:motion.idxe3);             
plot(xa,ya3)                                        % Yaw rate against time
title('Yaw rate')

subplot(2,2,4)
ya4 = flightdata.delta_a.data(motion.idx3:motion.idxe3);   % control input
plot(xa,ya4)
hold on
ya5 = flightdata.delta_r.data(motion.idx3:motion.idxe3);
plot(xa,ya5)
legend('delta a','delta r')
hold off
title('control input')

hold off

%Motion 4
motion.idx4 = find(flightdata.time.data==motion.t4);               % Starting time
motion.idxe4 = find(flightdata.time.data==motion.t4+45);              %temporary end time
refmotion(4) = figure(4);
title('Yaw damped Dutch roll')
subplot(2,2,1)
xa = flightdata.time.data(motion.idx4:motion.idxe4)-flightdata.time.data(motion.idx4);
ya1 = (pi/180).*flightdata.Ahrs1_Roll.data(motion.idx4:motion.idxe4);
plot(xa,ya1)                                        % Roll against time
title('Roll')

subplot(2,2,2)
ya2 = (pi/180).*flightdata.Ahrs1_bRollRate.data(motion.idx4:motion.idxe4);
plot(xa,ya2)                                        % Roll rate against time
title('Roll rate')

subplot(2,2,3)
ya3 = (pi/180).*flightdata.Ahrs1_bYawRate.data(motion.idx4:motion.idxe4);             
plot(xa,ya3)                                        % Yaw rate against time
title('Yaw rate')

subplot(2,2,4)
ya4 = flightdata.delta_a.data(motion.idx4:motion.idxe4);   % control input
plot(xa,ya4)
hold on
ya5 = flightdata.delta_r.data(motion.idx4:motion.idxe4);
plot(xa,ya5)
legend('delta a','delta r')
hold off
title('control input')

hold off

%Motion 5
motion.idx5 = find(flightdata.time.data==motion.t5);               % Starting time
motion.idxe5 = find(flightdata.time.data==motion.t5+50);              %temporary end time
refmotion(5) = figure(5);
title('Aperiodic Roll')
subplot(2,2,1)
xa = flightdata.time.data(motion.idx5:motion.idxe5)-flightdata.time.data(motion.idx5);
ya1 = (pi/180).*flightdata.Ahrs1_Roll.data(motion.idx5:motion.idxe5);
plot(xa,ya1)                                        % Roll against time
title('Roll')

subplot(2,2,2)
ya2 = (pi/180).*flightdata.Ahrs1_bRollRate.data(motion.idx5:motion.idxe5);
plot(xa,ya2)                                        % Roll rate against time
title('Roll rate')

subplot(2,2,3)
ya3 = (pi/180).*flightdata.Ahrs1_bYawRate.data(motion.idx5:motion.idxe5);             
plot(xa,ya3)                                        % Yaw against time
title('Yaw rate')

subplot(2,2,4)
ya4 = flightdata.delta_a.data(motion.idx5:motion.idxe5);   % control input
plot(xa,ya4)
hold on
ya5 = flightdata.delta_r.data(motion.idx5:motion.idxe5);
plot(xa,ya5)
legend('delta a','delta r')
hold off
title('control input')


hold off

%Motion 6
motion.idx6 = find(flightdata.time.data==motion.t6);               % Starting time
motion.idxe6 = find(flightdata.time.data==motion.t6+220);              %temporary end time
refmotion(6) = figure(6);
title('Spiral')
subplot(2,2,1)
xa = flightdata.time.data(motion.idx6:motion.idxe6)-flightdata.time.data(motion.idx6);
ya1 = (pi/180).*flightdata.Ahrs1_Roll.data(motion.idx6:motion.idxe6);
plot(xa,ya1)                                        % Roll against time
title('Roll')

subplot(2,2,2)
ya2 = (pi/180).*flightdata.Ahrs1_bRollRate.data(motion.idx6:motion.idxe6);
plot(xa,ya2)                                        % Roll rate against time
title('Roll rate')

subplot(2,2,3)
ya3 = (pi/180).*flightdata.Ahrs1_bYawRate.data(motion.idx6:motion.idxe6);            
plot(xa,ya3)                                        % Yaw rate against time
title('Yaw rate')

subplot(2,2,4)
ya4 = flightdata.delta_a.data(motion.idx6:motion.idxe6);   % control input
plot(xa,ya4)
hold on
ya5 = flightdata.delta_r.data(motion.idx6:motion.idxe6);
plot(xa,ya5)
legend('delta a','delta r')
hold off
title('control input')


hold off
savefig(refmotion,'referencemotions.fig')
% Phi = flightdata.Ahrs1_Roll.data;                   % Roll angle(deg)
% 
% p = flightdata.Ahrs1_bRollRate.data;                % Roll rate (deg/s)
% 
% %Psi = flightdata. It seems like this has to be calculated from flight heading    % Yaw(deg)
% 
% r = flightdata.Ahrs1_bYawRate.data;                 % Yaw rate(deg/s)


clearvars xa xs ya1 ya2 ya3 ya4 ya5 ys1 ys2 ys3 ys4 t1 t2 t3 t4 t5 idx1 idxe1 idx2 idxe2 idx3 idxe3 idx4 idxe4 idx5 idxe5
















