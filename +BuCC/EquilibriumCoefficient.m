classdef EquilibriumCoefficient < Geochemistry_Helpers.pX
    properties
        conditions;
        
        pressure_correction = [NaN,NaN,NaN,NaN,NaN]
        correction = NaN;
        scaled_correction = NaN;
        function_handle = NaN;
        function_coefficients = NaN;
        
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
        function self = EquilibriumCoefficient(name);
            self.conditions = BuCC.Conditions();
            if nargin>0
                try
                    bucc_search = what("+Bucc");
                    current_directory = pwd;
                    bucc_search = strrep(strrep(join([bucc_search(1).path,"\Configuration"],""),current_directory,"."),"\","/");
                    
                    
                    raw_pressure_file_contents = fileread(bucc_search+"/equilibrium_coefficient_pressure_correction.json");
                    json_pressure_file_contents = jsondecode(raw_pressure_file_contents);
                    valid_json_file_found = 1;
                    
                    self.pressure_correction = json_pressure_file_contents.(name);
                    
                    raw_function_file_contents = fileread(BuCC_search(1).path+"/Configuration/equilibrium_coefficient_functions.json");
                    json_function_file_contents = jsondecode(raw_function_file_contents);
                    self.function_handle = str2func(json_function_file_contents.(name));
                catch
                    valid_json_file_found = 0;
                    self.pressure_correction = [NaN,NaN,NaN,NaN,NaN];
                    self.function_handle = NaN;
                end
            end
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
        function output = get.temperature(self)
            output = self.conditions.temperature;
        end
        function output = get.salinity(self)
            output = self.conditions.salinity;
        end
        function output = get.pressure(self)
            output = self.conditions.pressure;
        end
        function output = get.calcium(self)
            output = self.conditions.calcium;
        end
        function output = get.magnesium(self)
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
        function calculate(self)
            ionic_strength = (19.924*self.conditions.salinity)/(1000-1.005*self.conditions.salinity);
            self.value = exp(self.function_handle(self.function_coefficients,self.conditions.temperature+273.15,self.conditions.salinity,ionic_strength));
        end
    end
end