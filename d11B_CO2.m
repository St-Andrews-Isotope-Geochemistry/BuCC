classdef d11B_CO2 < handle&Collator
    properties
        species_calibration
        boron
        carbonate_chemistry
    end
    methods
        % Constructor
        function self = d11B_CO2()
            % addpath("./Geochemistry_Helpers");
%             addpath("./Boron_Species_Calibration");
%             addpath("./Boron_Systematics");
%             addpath("./CO2_Systematics");
%             addpath("./CO2_Systematics/MyAMI");
%             addpath(genpath("./CO2_Systematics"));
            
            self.species_calibration = BoronSpeciesCalibration("polynomial",[1,0]);
            self.boron = Boron_pH();
            self.carbonate_chemistry = CarbonateChemistry();
        
            % Link together
            self.boron.d11B_4 = self.species_calibration.d11B_4;
            self.boron.pKb = self.carbonate_chemistry.equilibrium_coefficients.kb;
            self.boron.pH = self.carbonate_chemistry.pH; 
        end
        
        function calculate(self)
            for self_index = 1:numel(self);
                self(self_index).species_calibration.apply();
                self(self_index).carbonate_chemistry.equilibrium_coefficients.calculate();
                self(self_index).boron.calculate();
                self(self_index).carbonate_chemistry.calculate();
            end
        end
    end
end