classdef EquilibriumCoefficient < pX
    properties
        conditions;
        
        pressure_correction = [NaN,NaN,NaN,NaN,NaN]
        correction = NaN;
        scaled_correction = NaN;
        
        R_gas = 8.314510;        % J mol-1 deg-1 (perfect Gas)  
        R_P = 83.14472;             % mol bar deg-1 
    end
    properties (Dependent=true)
        temperature
        salinity
        pressure
        calcium
        magnesium
    end
    properties (Hidden=true)
        value_no_pressure_correction = NaN;
    end
    methods
        % Constructor
        function self = EquilibriumCoefficient();
            self.conditions = Conditions();
        end
        
        % Setters
        function set.temperature(self,value)
            self.conditions.temperature = value;
        end
        function set.salinity(self,value)
            self.conditions.salinity = value;
        end
        function set.pressure(self,value)
            self.conditions.pressure = value;
        end
        function set.calcium(self,value)
            self.conditions.calcium = value;
        end
        function set.magnesium(self,value)
            self.conditions.magnesium = value;
        end
        
        % Getters
        function output = get.temperature(self);
            output = self.conditions.temperature;
        end
        function output = get.salinity(self);
            output = self.conditions.salinity;
        end
        function output = get.pressure(self);
            output = self.conditions.pressure;
        end
        function output = get.calcium(self);
            output = self.conditions.calcium;
        end
        function output = get.magnesium(self);
            output = self.conditions.magnesium;
        end
        
        function doPressureCorrection(self)
            deltaV = self.pressure_correction(1) + self.pressure_correction(2)*self.temperature + self.pressure_correction(3)*self.temperature^2;
            deltaK = self.pressure_correction(4) + self.pressure_correction(5)*self.temperature;
            self.correction = exp(-(deltaV./(self.R_P*(self.temperature+273.15)))*self.pressure + (0.5*deltaK/(self.R_P*(self.temperature+273.15)))*self.pressure^2);
            self.value = self.value*self.correction;
        end
        function doScaleCorrectedPressureCorrection(self,scale_correction)
            deltaV = self.pressure_correction(1) + self.pressure_correction(2)*self.temperature + self.pressure_correction(3)*self.temperature^2;
            deltaK = self.pressure_correction(4) + self.pressure_correction(5)*self.temperature;
            self.correction = exp(-(deltaV./(self.R_P*(self.temperature+273.15)))*self.pressure + (0.5*deltaK/(self.R_P*(self.temperature+273.15)))*self.pressure^2);
            self.scaled_correction = self.correction*scale_correction;
            self.value = self.value*self.scaled_correction;
        end
    end
end