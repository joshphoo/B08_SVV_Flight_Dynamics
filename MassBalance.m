OEW = 9165.0*0.45359237; %kg
OEWarm = 292.18*0.0254; %m
blockfuel = 2750; % block fuel veranderen bij real data
%-------------------Reading Stationary Data--------------------------------

%-------------------Put The Data In array----------------------------------

xcgdat = [131, 131, 170, 214, 214, 251, 251, 288, 288]'; %inch
fuelmoments = [100,298.16;200,591.18;300,879.08;400,1165.42;500,1448.40;600,1732.53;700,2014.80;800,2298.84;900,2581.92;1000,2866.30;1100,3150.18;1200,3434.52;1300,3718.52;1400,4003.23;1500,4287;1600,4572.24;1700,4856.56;1800,5141.16;1900,5425.64;2000,5709.90;2100,5994.04;2200,6278.47;2300,6562.82;2400,6846.96;2500,7131.00;2600,7415.33;2700,7699.60;2800,7984.34;2900,8269.06;3000,8554.05;3100,8839.04;3200,9124.80;3300,9410.62;3400,9696.97;3500,9983.40;3600,10270.08;3700,10556.84;3800,10843.87;3900,11131.00;4000,11418.20;4100,11705.50;4200,11993.31;4300,12281.18;4400,12569.04;4500,12856.86;4600,13144.73];
%plot(fuelmoments(:,1),fuelmoments(:,2))

% for i = 1:9
%     xcgdatummom = OEW*OEWarm + excel_data_reader.weights(i)*(0.0254*xcgdat(i)); %moment [kg*m]
% end

xcgdatummom = OEW*OEWarm + sum(excel_data_reader.weights.*xcgdat*0.0254);

Fweight1 = blockfuel - excel_data_reader.Fused1; %lbs
Fweight2 = blockfuel - excel_data_reader.Fused2;
Fweight3 = blockfuel - excel_data_reader.Fused3;

%xcg1
for i = 1:length(Fweight1)
    for j = 1:46
        if (Fweight1(i) >=fuelmoments(j,1) ) && (Fweight1(i) <=(fuelmoments(j,1)+100))
            fuelmom1 = ((Fweight1-fuelmoments(j,1))*((fuelmoments((j+1),2)-fuelmoments(j,2))/100)+fuelmoments(j,2))*0.0254*0.45359237*100;
        end
    end
end 
xcg1 = []; %[kg*m]

for i = 1:length(Fweight1)
    massbalance.Weight1 = OEW + sum(excel_data_reader.weights) + Fweight1.*0.45359237; %kg
    xcg1 = [xcg1,-(((xcgdatummom + fuelmom1(i))/massbalance.Weight1(i))-6.643624)]; %m, - teken betekent ... meter achter de MAC zodat x positief naar voren wijst vanaf mac
end
massbalance.xcg1 = xcg1;
%xcg2
for i = 1:length(Fweight2)
    for j = 1:46
        if (Fweight2(i) >=fuelmoments(j,1) ) && (Fweight2(i) <=(fuelmoments(j,1)+100))
            fuelmom2 = ((Fweight2-fuelmoments(j,1))*((fuelmoments((j+1),2)-fuelmoments(j,2))/100)+fuelmoments(j,2))*0.0254*0.45359237*100;
        end
    end
end 
xcg2 = [];
for i = 1:length(Fweight2)
    massbalance.Weight2 = OEW + sum(excel_data_reader.weights) + Fweight2.*0.45359237; 
    xcg2 = [xcg2,-(((xcgdatummom + fuelmom2(i))/massbalance.Weight2(i))-6.643624)];
end
massbalance.xcg2 = xcg2;
%xcg3
for i = 1:length(Fweight3)
    for j = 1:46
        if (Fweight3(i) >=fuelmoments(j,1) ) && (Fweight3(i) <=(fuelmoments(j,1)+100))
            fuelmom3 = ((Fweight3-fuelmoments(j,1))*((fuelmoments((j+1),2)-fuelmoments(j,2))/100)+fuelmoments(j,2))*0.0254*0.45359237*100;
        end
    end
end 
massbalance.Weight3 = OEW + sum(excel_data_reader.weights) + Fweight3.*0.45359237;

xcg3 = [-(((xcgdatummom + fuelmom3(1))/massbalance.Weight3(1))-6.643624)];
massbalance.xcg3 = [xcg3, -((((xcgdatummom - excel_data_reader.weights(9)*(154*0.0254)) + fuelmom3(2))/massbalance.Weight3(2))-6.643624)];
clearvars filename Fmax fuelmom1 fuelmom2 fuelmom3 fuelmoments excel_data_reader.Fused1 ...
          excel_data_reader.Fused2 excel_data_reader.Fused3 Fweight1 Fweight2 Fweight3 i j OEW OEWarm ...
          excel_data_reader.weights xcg1 xcg2 xcg3 xcgdat ...
          xcgdatummom