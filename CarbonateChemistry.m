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
        
        conditions
        equilibrium_coefficients
        
        units = "xmol/kg";
        units_value = NaN;
        units_set = false;
    end
    properties (Dependent=true)
        temperature
        salinity
        pressure
        calcium
        magnesium
    end
    methods
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
        function set.units(self,value);
            self.units = value;
            self.getUnitsValue();
            self.units_set = true;
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
        
        % Constructor
        function self = CarbonateChemistry()
            self.conditions = Conditions();
            
            self.equilibrium_coefficients = EquilibriumCoefficients();
            self.equilibrium_coefficients.conditions = self.conditions;
            
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
            unit_normalisation = 10^self.units_value;
            
            if isnan(self.co2) && ~isnan(self.atmospheric_co2_partial_pressure)
                self.co2 = self.atmospheric_co2_partial_pressure*k0;
            elseif ~isnan(self.co2) && isnan(self.atmospheric_co2_partial_pressure)
                self.atmospheric_co2_partial_pressure = ((self.co2/unit_normalisation)/k0)*1e6;
            elseif isnan(self.co2) && isnan(self.atmospheric_co2_partial_pressure)
                atmospheric_co2_partial_pressure = self.co2/k0;
                if atmospheric_co2_partial_pressure~=self.atmospheric_co2_partial_pressure && ~isnan(atmospheric_co2_partial_pressure);
                    error("Inconsistent ocean CO2 and atmospheric CO2");
                end
            end
        end
        function calculate(self)
            for self_index=1:numel(self);                
                mgca_unit_normalisation = 10^self(self_index).conditions.mgca_units_value;
                calcium = self(self_index).calcium/mgca_unit_normalisation;
                magnesium = self(self_index).magnesium/mgca_unit_normalisation;
                
                if ~self(self_index).equilibrium_coefficients.calculated;
                    self(self_index).equilibrium_coefficients.calculate();
                end
                known_properties = self(self_index).getKnownProperties();

                if isnan(self(self_index).boron)
                    self(self_index).calculateBoronFromSalinity();
                end

                if ~isnan(self(self_index).co2) || ~isnan(self(self_index).atmospheric_co2_partial_pressure)
                    self(self_index).calculate_CO2();
                end
                if ~isnan(self(self_index).saturation_state)
                    self(self_index).co3 = ((self(self_index).saturation_state)*self(self_index).equilibrium_coefficients.kc.value)/(calcium);
                    known_properties(known_properties=="saturation_state") = "co3";
                end

                k1 = self(self_index).equilibrium_coefficients.k1.value;
                k2 = self(self_index).equilibrium_coefficients.k2.value;
                kw = self(self_index).equilibrium_coefficients.kw.value;
                kb = self(self_index).equilibrium_coefficients.kb.value;
            
                if ~isnan(self(self_index).pH.value)
                    pH = self(self_index).pH.value;
                    switch known_properties(2)
                        case "co2"
                            self(self_index).estimate_units("co2");
                            unit_normalisation = 10^self(self_index).units_value;

                            co2 = (self(self_index).co2/unit_normalisation);

                            dic = co2*(1+(k1/pH)+k1*(k2/pH^2));
                            hco3 = dic/(1+(pH/k1)+(k2/pH));
                            co3 = dic/(1+(pH/k2)+((pH^2)/(k1*k2)));
                            alkalinity = co2*((k1/pH)+((2*k1*k2)/pH^2))+(kb*self(self_index).boron)/(kb+pH)+(kw/pH)-pH;
                        case "hco3"
                            self(self_index).estimate_units("hco3");
                            unit_normalisation = 10^self(self_index).units_value;

                            hco3 = (self(self_index).hco3/unit_normalisation);

                            dic = hco3*(1+(pH/k1)+(k2/pH));
                            co2 = dic/(1+(k1/pH)+((k1*k2)/pH^2));
                            co3 = dic/(1+(pH/k2)+((pH^2)/(k1*k2)));
                            alkalinity = co2*((k1/pH)+((2*k1*k2)/pH^2))+(kb*self(self_index).boron)/(kb+pH)+(kw/pH)-pH;
                        case "co3"
                            self(self_index).estimate_units("co3");
                            unit_normalisation = 10^self(self_index).units_value;

                            co3 = (self(self_index).co3/unit_normalisation);

                            dic = co3*(1+(pH/k2)+(pH^2/(k1*k2)));
                            hco3 = dic/(1+(pH/k1)+(k2/pH));
                            co2 = dic/(1+(k1/pH)+((k1*k2)/pH^2));
                            alkalinity = co2*((k1/pH)+((2*k1*k2)/pH^2))+(kb*self(self_index).boron)/(kb+pH)+(kw/pH)-pH;
                        case "alkalinity"
                            self(self_index).estimate_units("alkalinity");
                            unit_normalisation = 10^self(self_index).units_value;

                            alkalinity = (self(self_index).alkalinity/unit_normalisation);

                            co2 = (alkalinity-((kb*self(self_index).boron)/(kb+pH))-(kw/pH)+pH)/((k1/pH)+2*((k1*k2)/pH^2));
                            dic = co2*(1+(k1/pH)+k1*(k2/pH^2));
                            hco3 = dic/(1+(pH/k1)+(k2/pH));
                            co3 = dic/(1+(pH/k2)+((pH^2)/(k1*k2)));
                        case "dic"
                            self(self_index).estimate_units("dic");
                            unit_normalisation = 10^self(self_index).units_value;

                            dic = (self(self_index).dic/unit_normalisation);

                            co2 = dic/(1+(k1/pH)+((k1*k2)/pH^2));
                            hco3 = dic/(1+(pH/k1)+(k2/pH));
                            co3 = dic/(1+(pH/k2)+((pH^2)/(k1*k2)));
                            alkalinity = co2*((k1/pH)+((2*k1*k2)/pH^2))+((kb*self(self_index).boron)/(kb+pH))+(kw/pH)-pH;
                        otherwise
                            error("Not implemented yet");
                    end
                    self(self_index).dic = dic*unit_normalisation;
                    self(self_index).alkalinity = alkalinity*unit_normalisation;
                    self(self_index).co2 = co2*unit_normalisation;
                    self(self_index).hco3 = hco3*unit_normalisation;
                    self(self_index).co3 = co3*unit_normalisation;

                    self(self_index).calculate_CO2();
                    self(self_index).saturation_state = (calcium*co3)/self(self_index).equilibrium_coefficients.kc.value;
                elseif ~isnan(self(self_index).dic) && ~isnan(self(self_index).co2)
                    self(self_index).estimate_units("dic");
                    unit_normalisation = 10^self(self_index).units_value;
                            
                    dic = (self(self_index).dic/unit_normalisation);
                    co2 = (self(self_index).co2/unit_normalisation);                    
                    
                    p2 = dic-co2;
                    p1 = -co2*k1;
                    p0 = -co2*k1*k2;
                    p = [p2,p1,p0];
                    r = roots(p);
                    pH = max(real(r));
                    
                    dic = co2*(1+(k1/pH)+k1*(k2/pH^2));
                    hco3 = dic/(1+(pH/k1)+(k2/pH));
                    co3 = dic/(1+(pH/k2)+((pH^2)/(k1*k2)));
                    alkalinity = co2*((k1/pH)+((2*k1*k2)/pH^2))+(kb*self(self_index).boron)/(kb+pH)+(kw/pH)-pH;
                        
                    self(self_index).dic = dic*unit_normalisation;
                    self(self_index).alkalinity = alkalinity*unit_normalisation;
                    self(self_index).co2 = co2*unit_normalisation;
                    self(self_index).hco3 = hco3*unit_normalisation;
                    self(self_index).co3 = co3*unit_normalisation;

%                     self(self_index).calculate_CO2();
                    self(self_index).saturation_state = ((calcium*co3)/self(self_index).equilibrium_coefficients.kc.value);
                    self(self_index).pH.value = pH;
                elseif ~isnan(self(self_index).alkalinity) && ~isnan(self(self_index).co2)
                    self(self_index).estimate_units("alkalinity");
                    unit_normalisation = 10^self(self_index).units_value;
                    
                    alkalinity = (self(self_index).alkalinity/unit_normalisation);
                    co2 = (self(self_index).co2/unit_normalisation);
                    
                    p4 = 1.;
                    p3 = kb+alkalinity;
                    p2 = alkalinity*kb-co2*k1-kb*self.boron-kw;
                    p1 = -co2*kb*k1-co2*2.*k1*k2-kw*kb;
                    p0 = -2.*co2*kb*k1*k2;
                    p = [p4,p3,p2,p1,p0];
                    r = roots(p);
                    pH = max(real(r));
                    
                    dic = co2*(1+(k1/pH)+k1*(k2/pH^2));
                    hco3 = dic/(1+(pH/k1)+(k2/pH));
                    co3 = dic/(1+(pH/k2)+((pH^2)/(k1*k2)));
                    
                    self(self_index).dic = dic*unit_normalisation;
                    self(self_index).alkalinity = alkalinity*unit_normalisation;
                    self(self_index).co2 = co2*unit_normalisation;
                    self(self_index).hco3 = hco3*unit_normalisation;
                    self(self_index).co3 = co3*unit_normalisation;
                    
                    %                     self(self_index).calculate_CO2();
                    self(self_index).saturation_state = ((calcium*co3)/self(self_index).equilibrium_coefficients.kc.value);
                    self(self_index).pH.value = pH;
                elseif ~isnan(self.co2) && ~isnan(self.co3)
                    self(self_index).estimate_units("co3");
                    co3_units = self(self_index).units_value;
                    self(self_index).units(1) = "x";
                    
                    self(self_index).estimate_units("co2");
                    co2_units = self(self_index).units_value;
                    assert(co3_units==co2_units,"Unit mismatch");
                    
                    unit_normalisation = 10^self(self_index).units_value;
                    
                    co3 = (self(self_index).co3/unit_normalisation);
                    co2 = (self(self_index).co2/unit_normalisation);                    
                    
                    p4 = -co3/k1/k2;
                    p3 = -co3/k2;
                    p2 = co2-co3;
                    p1 = co2*k1;
                    p0 = co2*k1*k2;
                    p = [p4,p3,p2,p1,p0];
                    r = roots(p);
                    pH = max(real(r));
                    
                    dic = co2*(1.+k1/pH+k1*k2/pH/pH);
                    hco3 = dic/(1+pH/k1+k2/pH);
                    alkalinity = co2*(k1/pH+2.*k1*k2/pH/pH)+kb*self.boron/(kb+pH)+kw/pH-pH;
                    
                    self(self_index).dic = dic*unit_normalisation;
                    self(self_index).alkalinity = alkalinity*unit_normalisation;
                    self(self_index).co2 = co2*unit_normalisation;
                    self(self_index).hco3 = hco3*unit_normalisation;
                    self(self_index).co3 = co3*unit_normalisation;
                    self(self_index).saturation_state = ((calcium*co3)/self(self_index).equilibrium_coefficients.kc.value);
                    self(self_index).pH.value = pH;
                else
                    error("Not implemented yet");
                end
            end
        end
    end
end