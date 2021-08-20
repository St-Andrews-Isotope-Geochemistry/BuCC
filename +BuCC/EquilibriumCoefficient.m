classdef EquilibriumCoefficient < Geochemistry_Helpers.pX
    properties
        name;
        conditions;
        MyAMI;
        
        pressure_correction = [NaN,NaN,NaN,NaN,NaN]
        correction = NaN;
        scaled_correction = NaN;
        function_handle = NaN;
        function_coefficients = NaN;
        
%         R_gas = 8.314510;        % J mol-1 deg-1 (perfect Gas)  
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
        valid_file_found = false;
    end
    methods
        % Constructor
        function self = EquilibriumCoefficient(name,container_flag)
            if nargin>0
                self.name = name;
            end
            if nargin<2
                container_flag = 0;
            end
            if ~container_flag
                self.conditions = BuCC.Conditions();
                try
                    if self.valid_file_found
                        json_pressure = jsondecode(fileread("equilibrium_coefficient_pressure_correction.json"));
                        json_function = jsondecode(fileread("equilibrium_coefficient_functions.json"));
                        
                        self.parsePressureCorrectionAndFunctions(json_pressure,json_function);
                        
                        self.valid_file_found = true;
                    end
                catch
                    valid_json_file_found = 0;
                    self.pressure_correction = [NaN,NaN,NaN,NaN,NaN];
                    self.function_handle = NaN;
                end
            end
        end
        function self = parsePressureCorrectionsAndFunctions(self,json_pressure,json_function)            
            self.pressure_correction = json_pressure.(self.name);
            self.function_handle = str2func(json_function.(self.name));
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
            output = self.conditions.oceanic_pressure;
        end
        function output = get.calcium(self)
            output = self.conditions.calcium;
        end
        function output = get.magnesium(self)
            output = self.conditions.magnesium;
        end
        
        function self = doPressureCorrection(self)
            deltaV = self.pressure_correction(1) + self.pressure_correction(2)*self.temperature + self.pressure_correction(3)*self.temperature^2;
            deltaK = self.pressure_correction(4) + self.pressure_correction(5)*self.temperature;
            self.correction = exp(-(deltaV./(self.R_P*(self.temperature+273.15)))*self.pressure + (0.5*deltaK/(self.R_P*(self.temperature+273.15)))*self.pressure^2);
            self.value = self.value*self.correction;
        end
        function self = doScaleCorrectedPressureCorrection(self,scale_correction)
            deltaV = self.pressure_correction(1) + self.pressure_correction(2)*self.temperature + self.pressure_correction(3)*self.temperature^2;
            deltaK = self.pressure_correction(4) + self.pressure_correction(5)*self.temperature;
            self.correction = exp(-(deltaV./(self.R_P*(self.temperature+273.15)))*self.pressure + (0.5*deltaK/(self.R_P*(self.temperature+273.15)))*self.pressure^2);
            self.scaled_correction = self.correction*scale_correction;
            self.value = self.value*self.scaled_correction;
        end
        function self = calculate(self)
            for self_index = 1:numel(self)
                if isempty(self(self_index).MyAMI)
                    self(self_index).MyAMI = MyAMI.MyAMI("Precalculated",true);
                end
                mgca_unit_normalisation = 10^self(self_index).conditions.mgca_units_value;
                self(self_index).MyAMI.calculate(self(self_index).conditions.temperature,self(self_index).conditions.salinity,self(self_index).conditions.calcium/mgca_unit_normalisation,self(self_index).conditions.magnesium/mgca_unit_normalisation);
%                 self(self_index).function_coefficients = self(self_index).MyAMI.results(self(self_index).name);

                self(self_index).value = self(self_index).MyAMI.results(self(self_index).name);
                self(self_index).doPressureCorrection();
            end
        end
    end
end