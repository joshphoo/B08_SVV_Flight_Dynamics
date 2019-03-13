%-------------------Engine part--------------------------------------------
%-------------------Bring Data to correct values---------------------------
%Calculations
%-------------------Dataset 1
T1 = excel_data_reader.T1;
T1(:,4) = T1(:,4).*0.3048;
%This converts ft to m

T1(:,5) = parameter.M1;
%obtains mach number from marloes parameters

T1(:,10) = T1(:,10) + 273.15;
%Convert Celcius to Kelvin
T1(:,10) = T1(:,10)-((273.15+15)+(T1(:,4).*(-0.0065)));

%Convert pounds per hour to kg/s
T1(:,7) = T1(:,7).*0.000125998;
T1(:,8) = T1(:,8).*0.000125998;


INPUT1 = T1(:,[4,5,10,7,8]);
%This formats the input for thurst.exe correct

%-------------------Dataset 2
T2 = excel_data_reader.T2;
T2(:,4) = T2(:,4).*0.3048;
%This converts ft to m

T2(:,5) = parameter.M2;
%obtains mach number from marloes parameters

T2(:,13) = T2(:,13) + 273.15;
%Convert Celcius to Kelvin
T2(:,13) = T2(:,13)-((273.15+15)+(T2(:,4).*(-0.0065)));

%Convert pounds per hour to kg/s
T2(:,10) = T2(:,10).*0.000125998;
T2(:,11) = T2(:,11).*0.000125998;


INPUT2 = T2(:,[4,5,13,10,11]);
%This formats the input for thurst.exe correct

%-------------------Dataset 3
T3 = excel_data_reader.T3;
T3(:,4) = T3(:,4).*0.3048;
%This converts ft to m

T3(:,5) = parameter.M3;
%obtains mach number from marloes parameters

T3(:,13) = T3(:,13) + 273.15;
%Convert Celcius to Kelvin
T3(:,13) = T3(:,13)-((273.15+15)+(T3(:,4).*(-0.0065)));

%Convert pounds per hour to kg/s
T3(:,10) = T3(:,10).*0.000125998;
T3(:,11) = T3(:,11).*0.000125998;


INPUT3 = T3(:,[4,5,13,10,11]);
%This formats the input for thurst.exe correct
%---------COMBINE ALL INPUT TABLES INTO ONE AND SAVE TO MATLAB.DAT---------
INPUT = [INPUT1;INPUT2;INPUT3];
dlmwrite('matlab.dat',INPUT,'delimiter',' ')
system('thrust.exe');
%gives thust for left and right engine

clearvars T1 T2 T3 INPUT1 INPUT2 INPUT3 INPUT