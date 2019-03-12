%-------------------EXCELL DATA READER-------------------------------------
%-------------------Reading Stationary Data--------------------------------
filename = fullfile('REFERENCE_Post_Flight_Datasheet_Flight.xlsx');
xlrange = 'A28:J33';
x2range = 'A59:M65';
x3range = 'A75:M76';
%-------------------Put The Data In array----------------------------------
%-------------------Dataset 1
[num1, text1, raw1] = xlsread(filename,xlrange);
num1(isnan(num1))=0;
text1 = str2double(text1);
text1(isnan(text1))=0;
text1 = [zeros(size(text1,1),1) text1];
format shortG
T1 = num1 + text1;

%-------------------Dataset 2
[num2, text2, raw2] = xlsread(filename,x2range);
num2(isnan(num2))=0;
num2 = [num2 zeros(size(num2,1),1)]
text2 = str2double(text2);
text2(isnan(text2))=0;
text2 = [zeros(size(text2,1),1) text2];
format shortG
num2
text2
T2 = num2 + text2;
%-------------------Dataset 3
[num3, text3, raw3] = xlsread(filename,x3range);
num3(isnan(num3))=0;
num3 = [num3 zeros(size(num3,1),1)]
text3 = str2double(text3);
text3(isnan(text3))=0;
text3 = [zeros(size(text3,1),1) text3];
format shortG
num3
text3
T3 = num3 + text3;