classdef d11BCO2 < handle&Geochemistry_Helpers.Collator
    properties
        species_calibration
        boron
        carbonate_chemistry
    end
    methods
        % Constructor
        function self = d11BCO2()            
            self.species_calibration = BuCC.BoronSpeciesCalibration("polynomial",[1,0]);
            self.boron = BuCC.BoronpH();
            self.carbonate_chemistry = BuCC.CarbonateChemistry();
        
            % Link together
            self.boron.d11B_4 = self.species_calibration.d11B_4;
            self.boron.pKb = self.carbonate_chemistry.equilibrium_coefficients.kb;
            self.boron.pH = self.carbonate_chemistry.pH; 
        end
        
        function calculate(self)
            for self_index = 1:numel(self)
                self(self_index).species_calibration.apply();
                self(self_index).carbonate_chemistry.equilibrium_coefficients.calculate();
                
                if self(self_index).boron.calculable
                    self(self_index).boron.calculate();
                    self(self_index).carbonate_chemistry.calculate();
                elseif self(self_index).carbonate_chemistry.calculable
                    self(self_index).carbonate_chemistry.calculate();
                    self(self_index).boron.calculate();   
                else
                    error("a");
                end
            end
        end
    end
    methods (Static=true)
        function output = create(number)
            output(number) = BuCC.d11BCO2();
            for index = 1:number
                output(index) = BuCC.d11BCO2();
            end
        end
    end
end