classdef d11B_CO2 < handle
    properties
        species_calibration
        boron
        carbonate_chemistry
    end
    methods
        % Constructor
        function self = d11B_CO2()
            addpath("./Boron_Species_Calibration");
            addpath("./Boron_Systematics");
            addpath("./CO2_Systematics");
%             addpath(genpath("./CO2_Systematics"));
            addpath("./Geochemistry_Helpers");
            
            self.species_calibration = BoronSpeciesCalibration();
            self.boron = Boron_pH();
            self.carbonate_chemistry = CarbonateChemistry();
        
            % Link together
            self.boron.d11B_4 = self.species_calibration.d11B_4;
            self.boron.pKb = self.carbonate_chemistry.equilibrium_coefficients.kb;
            self.boron.pH = self.carbonate_chemistry.pH; 
        end
        
        function calculate(self)
            self.species_calibration.apply();
            self.carbonate_chemistry.equilibrium_coefficients.calculate();
            self.boron.calculate();
            self.carbonate_chemistry.calculate();
        end
    end
end