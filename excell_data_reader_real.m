%-------------------EXCELL DATA READER-------------------------------------
%-------------------Read flightdata----------------------------------------
load('FTISxprt-20190306_143243.mat')
%-------------------Reading Stationary Data--------------------------------
filename = fullfile('Post_Flight_Datasheet_03_06_V4.xlsx');
x1range = 'A28:J33';
x1label = 'A25:J26';
x2range = 'A59:M64';
x2label = 'A56:M57';
x3range = 'A75:M76';
x3label = 'A72:M73';
%-------------------Put The Data In array----------------------------------
%-------------------Dataset 1
[num1, text1, raw1] = xlsread(filename,x1range);
[num,label1] = xlsread(filename,x1label);
num1(isnan(num1))=0;
num1 = [num1 zeros(size(text1,1),1)];
text1 = str2double(text1);
text1(isnan(text1))=0;
text1 = [zeros(size(text1,1),5) text1];
format shortG
T1 = num1 + text1;
excel_data_reader.T1 = T1;
excel_data_reader.label1 = label1;

%-------------------Dataset 2
[num2, text2, raw2] = xlsread(filename,x2range);
[num,label2] = xlsread(filename,x1label);
num2(isnan(num2))=0;
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
text3 = str2double(text3);
text3(isnan(text3))=0;
text3 = [zeros(size(text3,1),1) text3 zeros(size(text3,1),5)];
format shortG;
T3 = num3 + text3;
excel_data_reader.T3 = T3;
excel_data_reader.label3 = label3;


%------------------FOR MASS BALANCE

filename = fullfile('Post_Flight_Datasheet_03_06_V4.xlsx');
excel_data_reader.x1range1 = 'H8:H16';
excel_data_reader.weights = xlsread(filename,excel_data_reader.x1range1);
excel_data_reader.x2range2 = 'I28:I33';
excel_data_reader.Fused1 = xlsread(filename,excel_data_reader.x2range2);
excel_data_reader.x3range3 = 'L59:L64'; % naar flight data -> verander naar L64
excel_data_reader.Fused2 = xlsread(filename,excel_data_reader.x3range3);
excel_data_reader.x4range4 = 'L75:L76';
excel_data_reader.Fused3 = xlsread(filename,excel_data_reader.x4range4);

%------------------For graphs
motion.t1 = 0*3600 + 46*60 + 0;
motion.t2 = 0*3600 + 45*60 + 20;
motion.t3 = 0*3600 + 50*60 + 09;
motion.t4 = 0*3600 + 51*60 + 20;
motion.t5 = 0*3600 + 52*60 + 40;
motion.t6 = 0*3600 + 54*60 + 40;

motion.idx1 = find(flightdata.time.data==motion.t1-50);               % Starting time
motion.idxe1 = find(flightdata.time.data==motion.t1+310);              %temporary end time

motion.idx2 = find(flightdata.time.data==motion.t2-50);            % Starting time
motion.idxe2 = find(flightdata.time.data==motion.t2+150);              %temporary end time

motion.idx3 = find(flightdata.time.data==motion.t3-50);               % Starting time
motion.idxe3 = find(flightdata.time.data==motion.t3+119);              %temporary end time

motion.idx4 = find(flightdata.time.data==motion.t4-50);               % Starting time
motion.idxe4 = find(flightdata.time.data==motion.t4+145);              %temporary end time

motion.idx5 = find(flightdata.time.data==motion.t5-50);               % Starting time
motion.idxe5 = find(flightdata.time.data==motion.t5+150);              %temporary end time

motion.idx6 = find(flightdata.time.data==motion.t6-50);               % Starting time
motion.idxe6 = find(flightdata.time.data==motion.t6+260);              %temporary end time



clearvars filename x1range x2range x3range x1label x2label x3label num1 ...
        ` num2 num3 text1 text2 text3 raw1 raw2 raw3 T1 T2 T3 label1 ...
          label2 label3 num
