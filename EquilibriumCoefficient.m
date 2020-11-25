classdef EquilibriumCoefficient < pX
    properties
        temperature = NaN
        salinity = NaN
        pressure = NaN
        
        magnesium = NaN
        calcium = NaN
        
        pressure_correction = [NaN,NaN,NaN,NaN,NaN]
        correction = NaN;
        
        R_gas = 8.314510;        % J mol-1 deg-1 (perfect Gas)  
        R_P = 83.14472;             % mol bar deg-1 
    end
    properties (Hidden=true)
        value_no_pressure_correction = NaN;
    end
    methods
        function doPressureCorrection(self)
            deltaV = self.pressure_correction(1) + self.pressure_correction(2)*self.temperature + self.pressure_correction(3)*self.temperature^2;
            deltaK = self.pressure_correction(4) + self.pressure_correction(5)*self.temperature;
            self.correction = exp(-(deltaV./(self.R_P*(self.temperature+273.15)))*self.pressure + (0.5*deltaK/(self.R_P*(self.temperature+273.15)))*self.pressure^2);
            self.value = self.value*self.correction;
        end
    end
end