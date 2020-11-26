classdef BoronSpeciesCalibration < handle
    properties
        d11B_s = delta("B",NaN);
        d11B_4 = delta("B",NaN);
        
        form = "";
        coefficients = [];
    end
    methods
        function self = BoronSpeciesCalibration(form,coefficients)
            if nargin==0;
                self.form = "";
                self.coefficients = [];
            else
                self.form = form;
                self.coefficients = coefficients;
            end
        end
        function output = apply(self)
            if self.form=="polynomial" || self.form=="linear"
                output = 0;
                for index = 1:numel(self.coefficients)
                    output = output+self.coefficients(index)*self.d11B_s.value^(numel(self.coefficients)-index);
                end
                self.d11B_4.value = output;
            end
        end
        function assign(self,form,coefficients)
            self.form = form;
            self.coefficients = coefficients;
        end
    end
end