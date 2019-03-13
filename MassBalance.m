OEW = 9165.0*0.45359237;
OEWarm = 292.18*0.0254;
Fmax = 4050;
%-------------------Reading Stationary Data--------------------------------
filename = fullfile('REFERENCE_Post_Flight_Datasheet_Flight.xlsx');
xlrange1 = 'H8:H16';
weights = xlsread(filename,xlrange1);
x2range2 = 'I28:I33';
Fused = xlsread(filename,x2range2);
%-------------------Put The Data In array----------------------------------

xcgdat = [131, 131, 170, 214, 214, 251, 251, 288, 288];
fuelmoments = [100,298.16;200,591.18;300,879.08;400,1165.42;500,1448.40;600,1732.53;700,2014.80;800,2298.84;900,2581.92;1000,2866.30;1100,3150.18;1200,3434.52;1300,3718.52;1400,4003.23;1500,4287;1600,4572.24;1700,4856.56;1800,5141.16;1900,5425.64;2000,5709.90;2100,5994.04;2200,6278.47;2300,6562.82;2400,6846.96;2500,7131.00;2600,7415.33;2700,7699.60;2800,7984.34;2900,8269.06;3000,8554.05;3100,8839.04;3200,9124.80;3300,9410.62;3400,9696.97;3500,9983.40;3600,10270.08;3700,10556.84;3800,110843.87;3900,11131.00;4000,11418.20;4100,11705.50;4200,11993.31;4300,12281.18;4400,12569.04;4500,12856.86;4600,13144.73];
for i = 1:9
    xcgdatummom = OEW*OEWarm + weights(i)*(0.0254*xcgdat(i));
end
Fweight = 4050 - Fused;
fuelmom = [];
for i = 1:length(Fweight)
    for j = 1:46
        if (Fweight(i) >=fuelmoments(j,1) ) && (Fweight(i) <=(fuelmoments(j,1)+100))
            dudu = (Fweight-fuelmoments(j,1))*((fuelmoments((j+1),2)-fuelmoments(j,2))/100)+fuelmoments(j,2);
            fuelmom = [fuelmom,(0.011521246229174*dudu)];
        end
    end
end 
xcg = [];
for i = 1:length(Fweight)
    Weight = OEW + sum(weights) + (Fweight(i)*0.45359237);
    xcg = [xcg,-(((xcgdatummom + fuelmom(i))/Weight)-6.643624)];
end

        
