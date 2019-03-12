%-------------------Engine part--------------------------------------------
%-------------------Bring Data to correct values---------------------------
%Inputs:
gamma = 1.4;
%rho0 =
%p0 =
%Calculations
%-------------------Dataset 1
T1 = excel_data_reader.T1;
T1(:,4) = T1(:,4).*0.3048;
%This converts ft to m

%T(:,5) = T(:,5)-2
%This converts IAS to CAS

%T(;,5) = sqrt((2/gamma-1)*((1+(p0/p)((1+((gamma-1)*rho0*T(:,5).*T(:,5))/(2*gamma*p0))^(gamma/(gamma-1))-1))^((gamma-1)/gamma)-1))
%This converts CAS to M
T1(:,5) = 0.5*ones(6, 1);
%Temporarily sets Mach to 0.5

T1(:,10) = T1(:,10) + 273.15;
%Convert Celcius to Kelvin
T1(:,10) = T1(:,10)-((273.15+15)+(T1(:,4).*(-0.0065)));
disp('Check which temperature you have to subtract the Tisa from, total or static or a corrected thing?')
%Convert pounds per hour to kg/s
T1(:,7) = T1(:,7).*0.000125998;
T1(:,8) = T1(:,8).*0.000125998;


INPUT1 = T1(:,[4,5,10,7,8]);
%This formats the input for thurst.exe correct

%-------------------Dataset 2
T2 = excel_data_reader.T2;
T2(:,4) = T2(:,4).*0.3048;
%This converts ft to m

%T(:,5) = T(:,5)-2
%This converts IAS to CAS

%T(;,5) = sqrt((2/gamma-1)*((1+(p0/p)((1+((gamma-1)*rho0*T(:,5).*T(:,5))/(2*gamma*p0))^(gamma/(gamma-1))-1))^((gamma-1)/gamma)-1))
%This converts CAS to M
T2(:,5) = 0.5*ones((size(T2,1)), 1);
%Temporarily sets Mach to 0.5

T2(:,13) = T2(:,13) + 273.15;
%Convert Celcius to Kelvin
T2(:,13) = T2(:,13)-((273.15+15)+(T2(:,4).*(-0.0065)));
disp('Check which temperature you have to subtract the Tisa from, total or static or a corrected thing?')
%Convert pounds per hour to kg/s
T2(:,10) = T2(:,10).*0.000125998;
T2(:,11) = T2(:,11).*0.000125998;


INPUT2 = T2(:,[4,5,13,10,11]);
%This formats the input for thurst.exe correct

%-------------------Dataset 3
T3 = excel_data_reader.T3;
T3(:,4) = T3(:,4).*0.3048;
%This converts ft to m

%T(:,5) = T(:,5)-2
%This converts IAS to CAS

%T(;,5) = sqrt((2/gamma-1)*((1+(p0/p)((1+((gamma-1)*rho0*T(:,5).*T(:,5))/(2*gamma*p0))^(gamma/(gamma-1))-1))^((gamma-1)/gamma)-1))
%This converts CAS to M
T3(:,5) = 0.5*ones((size(T3,1)), 1);
%Temporarily sets Mach to 0.5

T3(:,13) = T3(:,13) + 273.15;
%Convert Celcius to Kelvin
T3(:,13) = T3(:,13)-((273.15+15)+(T3(:,4).*(-0.0065)));
disp('Check which temperature you have to subtract the Tisa from, total or static or a corrected thing?')
%Convert pounds per hour to kg/s
T3(:,10) = T3(:,10).*0.000125998;
T3(:,11) = T3(:,11).*0.000125998;


INPUT3 = T3(:,[4,5,13,10,11]);
%This formats the input for thurst.exe correct
%---------COMBINE ALL INPUT TABLES INTO ONE AND SAVE TO MATLAB.DAT---------
INPUT = [INPUT1;INPUT2;INPUT3];
dlmwrite('matlab.dat',INPUT,'delimiter',' ')
system('thrust.exe')
%gives thust for left and right engine
clearvars -except 'Input'


% % filename = str2double( regexp(YourCellArrayOfStrings, ',', '.') )
% % T1 = xlsread(filename,xlrange1);
% %writetable(T1,'matlab.dat','Delimiter',' ')
% [num, text, raw] = xlsread(filename,xlrange1)
% %T = str2num([num, text, raw])
%  [T1, cellarray] = xlsread(filename,xlrange1);        % Get Second Output
%  dbl_array = str2double(cellarray)                         
% 
% 
% %--------------------In case of needing to read flight data----------------
% % load('matlab.mat', 'flightdata')
% % T = table(flightdata.Dadc1_alt.data, flightdata.Dadc1_mach.data, flightdata.Dadc1_tat.data - flightdata.Dadc1_sat.data, flightdata.lh_engine_FMF.data, flightdata.rh_engine_FMF.data) ;
% % writetable(T,'matlab.dat','Delimiter',' ');
% % save matlab.dat
% % run thrust.exe