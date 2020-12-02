clear all

output_path = "./Results/";
output_name = "test_suite_object_oriented.mat";

temperature = 25;
salinity = 35;
pressure = 0;

calcium = 0.01;
magnesium = 0.05;

dic_index = 0;
pH_index = 0;
run_index = 0;

dic_options = [1500,2000,2500];
pH_options = [7.5,8.2,9.0];

run_index = run_index+1;

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


% At different temperature
dic_index = 0;
pH_index = 0;
run_index = run_index+1;
temperature = 10;
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

% At different salinity
dic_index = 0;
pH_index = 0;
run_index = run_index+1;
temperature = 25;
salinity = 40;
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

% At different pressure
dic_index = 0;
pH_index = 0;
run_index = run_index+1;
salinity = 35;
pressure = 10;
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

% At different calcium/magnesium
dic_index = 0;
pH_index = 0;
run_index = run_index+1;
pressure = 0;
calcium = 0.02;
magnesium = 0.03;
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

save(output_path+output_name,"alkalinity");