clear all

output_path = "./Results/";
output_name = "test_combo_object_oriented.mat";

temperature = 25;
salinity = 35;
pressure = 0;

calcium = 0.01;
magnesium = 0.05;

dic_index = 0;
pH_index = 0;
run_index = 0;
run_index = run_index+1;

dic_options = [1500,2000,2500];
pH_options = [7.5,8.2,9.0];

% pH + DIC
for dic = dic_options;
    dic_index = dic_index+1;
    pH_index = 0;
    for pH = pH_options;
        pH_index = pH_index+1;
        
        test = CarbonateChemistry();

        test.temperature = temperature;
        test.salinity = salinity;
        test.pressure = pressure;
        test.calcium = calcium;
        test.magnesium = magnesium;

        test.dic = dic_options(dic_index);
        test.pH.pValue = pH_options(pH_index);

        test.calculate();
        alkalinity(pH_index,dic_index,run_index) = test.alkalinity;
    end
end


% pH + co2
co2_index = 0;
pH_index = 0;
run_index = run_index+1;

co2_options = [5,10,20];
pH_options = [7.5,8.2,9.0];

for co2 = co2_options;
    co2_index = co2_index+1;
    pH_index = 0;
    for pH = pH_options;
        pH_index = pH_index+1;
        
        test = CarbonateChemistry();

        test.temperature = temperature;
        test.salinity = salinity;
        test.pressure = pressure;
        test.calcium = calcium;
        test.magnesium = magnesium;

        test.co2 = co2_options(co2_index);
        test.pH.pValue = pH_options(pH_index);

        test.calculate();
        alkalinity(pH_index,co2_index,run_index) = test.alkalinity;
    end
end

% pH + hco3
hco3_index = 0;
pH_index = 0;
run_index = run_index+1;

hco3_options = [1500,1800,2200];
pH_options = [7.5,8.2,9.0];

for hco3 = hco3_options;
    hco3_index = hco3_index+1;
    pH_index = 0;
    for pH = pH_options;
        pH_index = pH_index+1;
        
        test = CarbonateChemistry();

        test.temperature = temperature;
        test.salinity = salinity;
        test.pressure = pressure;
        test.calcium = calcium;
        test.magnesium = magnesium;

        test.hco3 = hco3_options(hco3_index);
        test.pH.pValue = pH_options(pH_index);

        test.calculate();
        alkalinity(pH_index,hco3_index,run_index) = test.alkalinity;
    end
end

% pH + co3
co3_index = 0;
pH_index = 0;
run_index = run_index+1;

co3_options = [50,100,200];
pH_options = [7.5,8.2,9.0];

for co3 = co3_options;
    co3_index = co3_index+1;
    pH_index = 0;
    for pH = pH_options;
        pH_index = pH_index+1;
        
        test = CarbonateChemistry();

        test.temperature = temperature;
        test.salinity = salinity;
        test.pressure = pressure;
        test.calcium = calcium;
        test.magnesium = magnesium;

        test.co3 = co3_options(co3_index);
        test.pH.pValue = pH_options(pH_index);

        test.calculate();
        alkalinity(pH_index,co3_index,run_index) = test.alkalinity;
    end
end

% pH + atmospheric co2
aco2_index = 0;
pH_index = 0;
run_index = run_index+1;

aco2_options = [50,100,200];
pH_options = [7.5,8.2,9.0];

for aco2 = aco2_options;
    aco2_index = aco2_index+1;
    pH_index = 0;
    for pH = pH_options;
        pH_index = pH_index+1;
        
        test = CarbonateChemistry();

        test.temperature = temperature;
        test.salinity = salinity;
        test.pressure = pressure;
        test.calcium = calcium;
        test.magnesium = magnesium;

        test.atmospheric_co2_partial_pressure = aco2_options(aco2_index);
        test.pH.pValue = pH_options(pH_index);

        test.calculate();
        alkalinity(pH_index,aco2_index,run_index) = test.alkalinity;
    end
end

save(output_path+output_name,"alkalinity");