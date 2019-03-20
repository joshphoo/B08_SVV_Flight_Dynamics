%-------------------EXCELL DATA READER-------------------------------------
%-------------------Reading Stationary Data--------------------------------
filename = fullfile('Post_Flight_Datasheet_03_06_V4.xlsx');
xlrange = 'A28:J33';
x1label = 'A25:J26';
x2range = 'A59:M64';
x2label = 'A56:M57';
x3range = 'A75:M76';
x3label = 'A72:M73';
%-------------------Put The Data In array----------------------------------
%-------------------Dataset 1
[num1, text1, raw1] = xlsread(filename,xlrange);
[num,label1] = xlsread(filename,x1label);
num1(isnan(num1))=0;
num1 = [num1 zeros(size(text1,1),1)];
text1 = str2double(text1);
text1(isnan(text1))=0;
text1 = [zeros(size(text1,1),5) text1];
format shortG
T1 = num1 + text1;
excel_data_reader_real.T1 = T1;
excel_data_reader_real.label1 = label1;

%-------------------Dataset 2
[num2, text2, raw2] = xlsread(filename,x2range);
[num,label2] = xlsread(filename,x1label);
num2(isnan(num2))=0;
text2 = str2double(text2);
text2(isnan(text2))=0;
text2 = [zeros(size(text2,1),1) text2];
format shortG;
T2 = num2 + text2;
excel_data_reader_real.T2 = T2;
excel_data_reader_real.label2 = label2;

%-------------------Dataset 3
[num3, text3, raw3] = xlsread(filename,x3range);
[num,label3] = xlsread(filename,x1label);
num3(isnan(num3))=0;
text3 = str2double(text3);
text3(isnan(text3))=0;
text3 = [zeros(size(text3,1),1) text3 zeros(size(text3,1),5)];
format shortG;
T3 = num3 + text3;
excel_data_reader_real.T3 = T3;
excel_data_reader_real.label3 = label3;

clearvars filename x1range x2range x3range x1label x2label x3label num1 ...
        ` num2 num3 text1 text2 text3 raw 1 raw2 raw3 T1 T2 T3