classdef BoronSpeciesCalibration < handle&Geochemistry_Helpers.Collator
    properties
        d11B_measured
        d11B_4
        
        form = "";
        coefficients = [];
    end
    properties (Dependent=true)
        calculable
    end
    methods
        function self = BoronSpeciesCalibration(form,coefficients)
            if nargin==0
                self.form = "";
                self.coefficients = [];
            else
                self.form = form;
                self.coefficients = coefficients;
            end
            self.d11B_measured = Geochemistry_Helpers.delta("B",NaN);
            self.d11B_4 = Geochemistry_Helpers.delta("B",NaN);
        end
        
        function output = get.calculable(self)
            if (~isnan(self.d11B_measured) || ~isnan(self.d11B_4)) && (~strcmp(self.form,"")) && ~isempty(self.coefficients)
                output = true;
            else
                output = false;
            end
        end
        
        function output = calculate(self)
            for self_index = 1:numel(self)
                if self(self_index).form=="polynomial" || self(self_index).form=="linear"
                    if isnan(self(self_index).d11B_4.value)
                        if self(self_index).form=="polynomial"
                            assert(numel(self(self_index).coefficients)==2,"Must be linear to be invertible");
                        end
                        output = (self(self_index).d11B_measured.value-self(self_index).coefficients(2))./self(self_index).coefficients(1);
                        self(self_index).d11B_4.value = output;
                    elseif isnan(self(self_index).d11B_measured.value)                        
                        output = 0;
                        for index = 1:numel(self(self_index).coefficients)
                            output = output+self(self_index).coefficients(index)*self(self_index).d11B_4.value^(numel(self(self_index).coefficients)-index);
                        end
                        self(self_index).d11B_measured.value = output;
                    else
                        error("Species calibration issue");
                    end
                end
            end
        end
        function assign(self,form,coefficients)
            self.form = form;
            self.coefficients = coefficients;
        end
    end
end