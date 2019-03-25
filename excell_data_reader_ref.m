%-------------------EXCELL DATA READER-------------------------------------
%-------------------Read flightdata----------------------------------------
load('matlab.mat')
%-------------------Reading Stationary Data--------------------------------
filename = fullfile('REFERENCE_Post_Flight_Datasheet_Flight.xlsx');
x1range = 'A28:J33';
x1label = 'A25:J26';
x2range = 'A59:M65';
x2label = 'A56:M57';
x3range = 'A75:M76';
x3label = 'A72:M73';
%-------------------Put The Data In array----------------------------------
%-------------------Dataset 1
[num1, text1, raw1] = xlsread(filename,x1range);
[num,label1] = xlsread(filename,x1label);
num1(isnan(num1))=0;
text1 = str2double(text1);
text1(isnan(text1))=0;
text1 = [zeros(size(text1,1),1) text1];
format shortG
T1 = num1 + text1;
excel_data_reader.T1 = T1;
excel_data_reader.label1 = label1;

%-------------------Dataset 2
[num2, text2, raw2] = xlsread(filename,x2range);
[num,label2] = xlsread(filename,x1label);
num2(isnan(num2))=0;
num2 = [num2 zeros(size(num2,1),1)];
text2 = str2double(text2);
text2(isnan(text2))=0;
text2 = [zeros(size(text2,1),1) text2];
format shortG;
T2 = num2 + text2;
excel_data_reader.T2 = T2;
excel_data_reader.label2 = label2;

%-------------------Dataset 3
[num3, text3, raw3] = xlsread(filename,x3range);
[num,label3] = xlsread(filename,x1label);
num3(isnan(num3))=0;
num3 = [num3 zeros(size(num3,1),1)];
text3 = str2double(text3);
text3(isnan(text3))=0;
text3 = [zeros(size(text3,1),1) text3];
format shortG;
T3 = num3 + text3;
excel_data_reader.T3 = T3;
excel_data_reader.label3 = label3;

clearvars filename x1range x2range x3range x1label x2label x3label num1 ...
        ` num2 num3 text1 text2 text3 raw1 raw2 raw3 T1 T2 T3 num ...
          label1 label2 label3
      
%-------------------FOR MASS BALANCE
filename = fullfile('REFERENCE_Post_Flight_Datasheet_Flight.xlsx');
excel_data_reader.x1range1 = 'H8:H16';
excel_data_reader.weights = xlsread(filename,excel_data_reader.x1range1);
excel_data_reader.x2range2 = 'I28:I33';
excel_data_reader.Fused1 = xlsread(filename,excel_data_reader.x2range2);
excel_data_reader.x3range3 = 'L59:L65'; % naar flight data -> verander naar L64
excel_data_reader.Fused2 = xlsread(filename,excel_data_reader.x3range3);
excel_data_reader.x4range4 = 'L75:L76';
excel_data_reader.Fused3 = xlsread(filename,excel_data_reader.x4range4);

%------------------For graphs
motion.t1 = 0*3600 + 53*60 + 56;
motion.t2 = 1*3600 + 0*60 + 33;
motion.t3 = 1*3600 + 1*60 + 57;
motion.t4 = 1*3600 + 2*60 + 47;
motion.t5 = 0*3600 + 59*60 + 10;
motion.t6 = 1*3600 + 5*60 + 38;

motion.idx1 = find(flightdata.time.data==motion.t1-10);               % Starting time
motion.idxe1 = find(flightdata.time.data==motion.t1+210);              %temporary end time

motion.idx2 = find(flightdata.time.data==motion.t2);            % Starting time
motion.idxe2 = find(flightdata.time.data==motion.t2+50);              %temporary end time

motion.idx3 = find(flightdata.time.data==motion.t3);               % Starting time
motion.idxe3 = find(flightdata.time.data==motion.t3+19);              %temporary end time

motion.idx4 = find(flightdata.time.data==motion.t4);               % Starting time
motion.idxe4 = find(flightdata.time.data==motion.t4+45);              %temporary end time

motion.idx5 = find(flightdata.time.data==motion.t5);               % Starting time
motion.idxe5 = find(flightdata.time.data==motion.t5+50);              %temporary end time

motion.idx6 = find(flightdata.time.data==motion.t6);               % Starting time
motion.idxe6 = find(flightdata.time.data==motion.t6+160);              %temporary end time
