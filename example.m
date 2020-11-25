clear

test = d11B_CO2();

test.species_calibration.assign("polynomial",[1,0]);
test.species_calibration.d11B_s.value = 20;

test.carbonate_chemistry.temperature = 25;
test.carbonate_chemistry.salinity = 25;
test.carbonate_chemistry.pressure = 0;
test.carbonate_chemistry.calcium = 0.01;
test.carbonate_chemistry.magnesium = 0.05;

test.carbonate_chemistry.dic = 2000;

test.calculate();