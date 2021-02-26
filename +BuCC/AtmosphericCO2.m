classdef AtmosphericCO2 < handle&Geochemistry_Helpers.Collator
    properties
        fugacity = NaN;
        partial_pressure = NaN;
        mole_fraction = NaN;
        conditions
        units = "x";
    end
    properties (Dependent=true,Hidden=true)
        f
        p
        x        
    end
    properties (Hidden=true)
        fugacity_units = " atm/atm";
        partial_pressure_units = " atm/atm";
        mole_fraction_units = "pp ";
        
        units_value
        units_set
    end
    methods
        % Constructor
        function self = AtmosphericCO2()
            self.conditions = BuCC.Conditions();
        end
        % Setters
        function set.f(self,value)
            self.fugacity = value;
        end
        function set.p(self,value)
            self.partial_pressure = value;
        end
        function set.x(self,value)
            self.mole_fraction = value;
        end
        function set.units(self,value)
            self.units = value;
            self.getUnitsValue();
            self.units_set = true;
        end
        
        % Getters
        function output = get.f(self)
            output = self.fugacity;
        end
        function output = get.p(self)
            output = self.partial_pressure;
        end
        function output = get.x(self)
            output = self.mole_fraction;
        end
        
        function output = isnan(self)
            if ~isnan(self.fugacity) || ~isnan(self.partial_pressure) || ~isnan(self.mole_fraction)
                output = false;
            else
                output = true;
            end
        end
        
        function estimateUnits(self)
            units_char = char(self.units);
            fugacity_units_char = char(self.fugacity_units);
            partial_pressure_units_char = char(self.partial_pressure_units);
            mole_fraction_units_char = char(self.mole_fraction_units);
            
            typical_values = 500e-6;
            
            % Use dic as the indicator case
            typical_order_of_magnitude = 3*floor(log10(typical_values)/3);
            input_order_of_magnitude = 3*floor(log10(self.mole_fraction)/3);
            difference_order_of_magnitude = input_order_of_magnitude-typical_order_of_magnitude;
            
            if units_char(1)=="x"
                self.units_value = difference_order_of_magnitude;
                
                if difference_order_of_magnitude==9
                    units_char(1) = "n";
                    fugacity_units_char(1) = "n";
                    partial_pressure_units_char(1) = "n";
                    mole_fraction_units_char(3) = "t";
                elseif difference_order_of_magnitude==6
                    units_char(1) = "μ";
                    fugacity_units_char(1) = "μ";
                    partial_pressure_units_char(1) = "μ";
                    mole_fraction_units_char(3) = "m";
                elseif difference_order_of_magnitude==3
                    units_char(1) = "m";                    
                    fugacity_units_char(1) = "m";
                    partial_pressure_units_char(1) = "m";
                    mole_fraction_units_char(3) = "k";
                elseif difference_order_of_magnitude==0
                    units_char(1) = " ";
                    fugacity_units_char(1) = " ";
                    partial_pressure_units_char(1) = " ";
                    mole_fraction_units_char(3) = "o";
                end
                self.units = units_char;
                self.fugacity_units = fugacity_units_char;
                self.partial_pressure_units = partial_pressure_units_char;
                self.mole_fraction_units = mole_fraction_units_char;
            end
        end
        function getUnitsValue(self)
            asChar = char(self.units);
            if asChar(1)=="m"
                self.units_value = 3;
            elseif asChar(1)=="u" || asChar(1)=="μ"
                self.units_value = 6;
            elseif asChar(1)==" "
                self.units_value = 0;
            else
                self.units_value = NaN;
            end            
        end
        
        function calculate(self)
            if ~isnan(self)
                temperature = self.conditions.temperature + 273.15;
                salinity = self.conditions.salinity;
                pressure = self.conditions.atmospheric_pressure;
                gas_constant = self.conditions.gas_constant;
                
                if ~isnan(self.partial_pressure)
                    self.mole_fraction = self.partialPressureToMoleFraction(self.partial_pressure,temperature,salinity,pressure);
                    self.fugacity = self.partialPressureToFugacity(self.partial_pressure,temperature,pressure,gas_constant);
                elseif ~isnan(self.mole_fraction)
                    self.partial_pressure = self.moleFractionToPartialPressure(self.mole_fraction,temperature,salinity,pressure);
                    self.fugacity = self.partialPressureToFugacity(self.partial_pressure,temperature,pressure,gas_constant);
                elseif ~isnan(self.fugacity)
                    self.partialPressure = self.fugacityToPartialPressure(self.fugacity,temperature,pressure,gas_constant);
                    self.mole_fraction = self.partialPressureToMoleFraction(self.partial_pressure,temperature,salinity,pressure);
                end
                
                self.estimateUnits();
            end
        end
        
        % Display
        function show(self,parameter)
            if strcmp(parameter,"partial_pressure") || strcmp(parameter,"p")
                disp(num2str(self.(parameter))+" "+string(self.partial_pressure_units));
            elseif strcmp(parameter,"mole_fraction") || strcmp(parameter,"x")
                disp(num2str(self.(parameter))+" "+string(self.mole_fraction_units));
            elseif strcmp(parameter,"fugacity") || strcmp(parameter,"f")
                disp(num2str(self.(parameter))+" "+string(self.fugacity_units));
            else
                error("Unsupported parameter specified");
            end
        end
    end
    methods (Static=true)
        function fugacity = partialPressureToFugacity(partial_pressure,temperature,pressure,gas_constant)
            co2_virial_coefficient = (-1636.75 + 12.0408*temperature - 3.27957e-2 * temperature^2 + 3.16528e-5 * temperature^3);
            delta = 57.7-0.0118*temperature;
            
            fugacity = partial_pressure*exp((101325/1e6)*pressure*((co2_virial_coefficient+2*delta)/(gas_constant*temperature)));
        end
        function partial_pressure = fugacityToPartialPressure(fugacity,temperature,pressure,gas_constant)
            co2_virial_coefficient = (-1636.75 + 12.0408*temperature - 3.27957e-2 * temperature^2 + 3.16528e-5 * temperature^3);
            delta = 57.7-0.0118*temperature;
            
            partial_pressure = fugacity/exp((101325/1e6)*pressure*((co2_virial_coefficient+2*delta)/(gas_constant*temperature)));
        end
        function partial_pressure = moleFractionToPartialPressure(mole_fraction,temperature,salinity,pressure)
            pH2O = exp(24.4543 - 6745.09/temperature - 4.8489*log(temperature/100) - 0.000544*salinity);
            partial_pressure = mole_fraction*(pressure-pH2O);
        end
        function mole_fraction = partialPressureToMoleFraction(partial_pressure,temperature,salinity,pressure)
            pH2O = exp(24.4543 - 6745.09/temperature - 4.8489*log(temperature/100) - 0.000544*salinity);
            mole_fraction = partial_pressure/(pressure-pH2O);
        end
    end
end