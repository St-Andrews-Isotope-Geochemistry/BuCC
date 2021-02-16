classdef Conditions < handle&Geochemistry_Helpers.Collator
    properties
        temperature = NaN;
        salinity = NaN;
        pressure = NaN;
        calcium = NaN;
        magnesium = NaN;
        
        
        mgca_units = "xmol/kg";
        mgca_units_value = NaN;
        mgca_units_set = false;
    end
    methods
        function self = Conditions(varargin)
            parser = inputParser;
            parser.KeepUnmatched = true;
            properties = ["temperature","salinity","pressure","calcium","magnesium"];
            
            for property = properties
                addOptional(parser,property,NaN);
            end
            
            parse(parser,varargin{:});
            
            for property = properties
                if ~isempty(parser.Results.(property)) && ~isnan(parser.Results.(property))
                    self.(property) = parser.Results.(property);
                end
            end
        end
        
        function set.calcium(self,value)
            self.calcium = value;
            if ~self.mgca_units_set
                self.estimate_units("calcium");
            end
        end
        function set.magnesium(self,value)
            self.magnesium = value;
            self.estimate_units("magnesium");
        end
        function set.mgca_units(self,value)
            self.mgca_units = value;
            self.getUnitsValue();
            self.mgca_units_set = true;
        end
        function estimate_units(self,parameter)
            units_char = char(self.mgca_units);
            
            typical_values_map = containers.Map(["calcium","magnesium"],[0.01,0.05]);
            
            % Use dic as the indicator case
            typical_order_of_magnitude = 3*floor(log10(typical_values_map(parameter))/3);
            input_order_of_magnitude = 3*floor(log10(self.(parameter))/3);
            difference_order_of_magnitude = input_order_of_magnitude-typical_order_of_magnitude;
            
            if units_char(1)=="x";
                self.mgca_units_value = difference_order_of_magnitude;
                
                if difference_order_of_magnitude==6;
                    units_char(1) = "μ";
                elseif difference_order_of_magnitude==3;
                    units_char(1) = "m";
                elseif difference_order_of_magnitude==0;
                    units_char(1) = " ";
                end
                self.mgca_units = units_char;
            elseif self.mgca_units_value ~= difference_order_of_magnitude;
                error("Different units for Ca and Mg not allowed");
            end
        end
        function getUnitsValue(self)
            asChar = char(self.mgca_units);
            if asChar(1)=="m"
                self.mgca_units_value = 3;
            elseif asChar(1)=="u" || asChar(1)=="μ"
                self.mgca_units_value = 6;
            elseif asChar(1)==" "
                self.mgca_units_value = 0;
            else
                self.mgca_units_value = NaN;
            end            
        end
    end
end