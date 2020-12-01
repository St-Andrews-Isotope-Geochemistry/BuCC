classdef CarbonateChemistry < handle&Collator
    properties
        pH
        
        dic = NaN;
        alkalinity = NaN;
        
        co2 = NaN;
        hco3 = NaN;
        co3 = NaN;
        atmospheric_co2_partial_pressure = NaN;
        
        boron = NaN
        
        saturation_state = NaN;
        
        temperature
        salinity
        pressure
        calcium
        magnesium
        
        equilibrium_coefficients
        
        units = "xmol/kg";
        units_value = NaN;
        units_set = false;
    end
    methods
        % Setters
        function set.temperature(self,value)
            self.temperature = value;
            self.equilibrium_coefficients.temperature = value;
        end
        function set.salinity(self,value)
            self.salinity = value;
            self.equilibrium_coefficients.salinity = value;
        end
        function set.pressure(self,value)
            self.pressure = value;
            self.equilibrium_coefficients.pressure = value;
        end        
        function set.calcium(self,value)
            self.calcium = value;
            self.equilibrium_coefficients.calcium = value;
        end        
        function set.magnesium(self,value)
            self.magnesium = value;
            self.equilibrium_coefficients.magnesium = value;
        end
        function set.units(self,value);
            self.units = value;
            self.getUnitsValue();
            self.units_set = true;
        end
        
        % Constructor
        function self = CarbonateChemistry()
            self.equilibrium_coefficients = EquilibriumCoefficients();
            self.pH = pX();
        end
        
        function calculateBoronFromSalinity(self)
            self.boron = (0.0002414/10.811)*(self.salinity/1.80655);
        end
        function known_properties = getKnownProperties(self)
            known_properties = [];
            
            if ~isnan(self.pH.value);
                known_properties = [known_properties,"pH"];
            end
            if ~isnan(self.dic);
                known_properties = [known_properties,"dic"];
            end
            if ~isnan(self.alkalinity);
                known_properties = [known_properties,"alkalinity"];
            end
            if ~isnan(self.co2) || ~isnan(self.atmospheric_co2_partial_pressure);
                known_properties = [known_properties,"co2"];
            end
            if ~isnan(self.hco3);
                known_properties = [known_properties,"hco3"];
            end
            if ~isnan(self.co3);
                known_properties = [known_properties,"co3"];
            end
            if ~isnan(self.saturation_state);
                known_properties = [known_properties,"saturation_state"];
            end            
        end
        
        function estimate_units(self,parameter);
            units_char = char(self.units);
            if units_char(1)=="x"
                typical_values_map = containers.Map(["dic","alkalinity","co2","hco3","co3"],[2000e-6,2300e-6,10e-6,1800e-6,200e-6]);

                % Use dic as the indicator case
                typical_order_of_magnitude = 3*floor(log10(typical_values_map(parameter))/3);
                input_order_of_magnitude = 3*floor(log10(self.(parameter))/3);
                difference_order_of_magnitude = input_order_of_magnitude-typical_order_of_magnitude;
                self.units_value = difference_order_of_magnitude;

                if difference_order_of_magnitude==6;
                    units_char(1) = "μ";
                elseif difference_order_of_magnitude==3;
                    units_char(1) = "m";
                elseif difference_order_of_magnitude==0;
                    units_char(1) = " ";
                end
                self.units = units_char;
            end
        end
        function getUnitsValue(self);
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
        
        function show(self,parameter);
            if parameter=="alkalinity" || parameter=="dic" || parameter=="co2" || parameter=="hco3" || parameter=="co3"
                disp(join([num2str(self.(parameter)),self.units]));
            else
                disp(self.(parameter));
            end
        end
        
        function calculate_CO2(self)
            k0 = self.equilibrium_coefficients.k0.value;
            
            if isnan(self.co2) && ~isnan(self.atmospheric_co2_partial_pressure)
                self.co2 = self.atmospheric_co2_partial_pressure*k0;
            elseif ~isnan(self.co2) && isnan(self.atmospheric_co2_partial_pressure)
                self.atmospheric_co2_partial_pressure = self.co2/k0;
            elseif isnan(self.co2) && isnan(self.atmospheric_co2_partial_pressure)
                atmospheric_co2_partial_pressure = self.co2/k0;
                if atmospheric_co2_partial_pressure~=self.atmospheric_co2_partial_pressure && ~isnan(atmospheric_co2_partial_pressure);
                    error("Inconsistent ocean CO2 and atmospheric CO2");
                end
            end
        end
        function calculate(self)
            if ~self.equilibrium_coefficients.calculated;
                self.equilibrium_coefficients.calculate();
            end
            known_properties = self.getKnownProperties();
            
            if isnan(self.boron)
                self.calculateBoronFromSalinity();
            end
            
            if ~isnan(self.co2) || ~isnan(self.atmospheric_co2_partial_pressure)
                self.calculate_CO2();
            end            
                
            k1 = self.equilibrium_coefficients.k1.value;
            k2 = self.equilibrium_coefficients.k2.value;
            kw = self.equilibrium_coefficients.kw.value;
            kb = self.equilibrium_coefficients.kb.value;
            
            if ~isnan(self.pH.value)
                pH = self.pH.value;
                switch known_properties(2)
                    case "co2"
                        self.estimate_units("co2");
                        unit_normalisation = 10^self.units_value;
                        
                        co2 = (self.co2/unit_normalisation);
                        
                        dic = co2*(1+(k1/pH)+k1*(k2/pH^2));
                        hco3 = dic/(1+(pH/k1)+(k2/pH));
                        co3 = dic/(1+(pH/k2)+((pH^2)/(k1*k2)));
                        alkalinity = co2*((k1/pH)+((2*k1*k2)/pH^2))+(kb*self.boron)/(kb+pH)+(kw/pH)-pH;
                    case "hco3"
                        self.estimate_units("hco3");
                        unit_normalisation = 10^self.units_value;
                        
                        hco3 = (self.hco3/unit_normalisation);
                        
                        dic = hco3*(1+(pH/k1)+(k2/pH));
                        co2 = dic/(1+(k1/pH)+((k1*k2)/pH^2));
                        co3 = dic/(1+(pH/k2)+((pH^2)/(k1*k2)));
                        alkalinity = co2*((k1/pH)+((2*k1*k2)/pH^2))+(kb*self.boron)/(kb+pH)+(kw/pH)-pH;
                    case "co3"
                        self.estimate_units("co3");
                        unit_normalisation = 10^self.units_value;
                        
                        co3 = (self.co3/unit_normalisation);
                        
                        dic = co3*(1+(pH/k2)+(pH^2/(k1*k2)));
                        hco3 = dic/(1+(pH/k1)+(k2/pH));
                        co2 = dic/(1+(k1/pH)+((k1*k2)/pH^2));
                        alkalinity = co2*((k1/pH)+((2*k1*k2)/pH^2))+(kb*self.boron)/(kb+pH)+(kw/pH)-pH;
                    case "alkalinity"
                        self.estimate_units("alkalinity");
                        unit_normalisation = 10^self.units_value;
                        
                        alkalinity = (self.alkalinity/unit_normalisation);
                        
                        co2 = (alkalinity-((kb*self.boron)/(kb+pH))-(kw/pH)+pH)/((k1/pH)+2*((k1*k2)/pH^2));
                        dic = co2*(1+(k1/pH)+k1*(k2/pH^2));
                        hco3 = dic/(1+(pH/k1)+(k2/pH));
                        co3 = dic/(1+(pH/k2)+((pH^2)/(k1*k2)));
                    case "dic"
                        self.estimate_units("dic");
                        unit_normalisation = 10^self.units_value;
                        
                        dic = (self.dic/unit_normalisation);
                        
                        co2 = dic/(1+(k1/pH)+((k1*k2)/pH^2));
                        hco3 = dic/(1+(pH/k1)+(k2/pH));
                        co3 = dic/(1+(pH/k2)+((pH^2)/(k1*k2)));
                        alkalinity = co2*((k1/pH)+((2*k1*k2)/pH^2))+(kb*self.boron)/(kb+pH)+(kw/pH)-pH;
                end
                self.dic = dic*unit_normalisation;
                self.alkalinity = alkalinity*unit_normalisation;
                self.co2 = co2*unit_normalisation;
                self.hco3 = hco3*unit_normalisation;
                self.co3 = co3*unit_normalisation;
                
                self.calculate_CO2();
                self.saturation_state = ((self.calcium*co3)/self.equilibrium_coefficients.kc.value);
            else
                error("Not implemented yet");
            end
        end
    end
end