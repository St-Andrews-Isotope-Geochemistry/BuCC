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
        
        % Constructor
        function self = CarbonateChemistry()
            self.equilibrium_coefficients = EquilibriumCoefficients();
            self.pH = pX();
        end
        
        function calculateBoronFromSalinity(self)
            self.boron = (0.0002414/10.811)*(self.salinity/1.80655);
        end
        function calculate_CO2(self)
            k0 = self.equilibrium_coefficients.k0.value;
            
            if isnan(self.co2) && ~isnan(self.atmospheric_co2_partial_pressure)
                self.co2 = self.atmospheric_co2_partial_pressure*k0;
            elseif ~isnan(self.co2) && isnan(self.atmospheric_co2_partial_pressure)
                self.atmospheric_co2_partial_pressure = self.co2/k0;
            elseif isnan(self.co2) && isnan(self.atmospheric_co2_partial_pressure)
                atmospheric_co2_partial_pressure = self.co2/k0;
                if atmospheric_co2_partial_pressure~=self.atmospheric_co2_partial_pressure
                    error("Inconsistent ocean CO2 and atmospheric CO2");
                end
            end
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
                switch known_properties(2)
                    case "co2"
                        self.dic = self.co2*(1+(k1/self.pH.value)+k1*(k2/self.pH.value^2));
                        self.hco3 = self.dic/(1+(self.pH.value/k1)+(k2/self.pH.value));
                        self.co3 = self.dic/(1+(self.pH.value/k2)+((self.pH.value^2)/(k1*k2)));
                        self.alkalinity = self.co2*((k1/self.pH.value)+((2.*k1*k2)/self.pH.value^2))+(kb*self.boron)/(kb+self.pH.value)+(kw/self.pH.value)-self.pH.value;
                    case "hco3"
                        self.dic = self.hco3*(1+(self.pH.value/k1)+(k2/self.pH.value));
                        self.co2 = self.dic/(1+(k1/self.pH.value)+((k1*k2)/self.pH.value.^2));
                        self.co3 = self.dic/(1+(self.pH.value/k2)+((self.pH.value^2)/(k1*k2)));
                        self.alkalinity = self.co2*((k1/self.pH.value)+((2.*k1*k2)/self.pH.value^2))+(kb*self.boron)/(kb+self.pH.value)+(kw/self.pH.value)-self.pH.value;
                    case "co3"
                        self.dic = self.co3*(1+(self.pH.value/k2)+(self.pH.value.^2/(k1*k2)));
                        self.hco3 = self.dic/(1+(self.pH.value/k1)+(k2/self.pH.value));
                        self.co2 = self.dic/(1+(k1/self.pH.value)+((k1*k2)/self.pH.value.^2));
                        self.alkalinity = self.co2*((k1/self.pH.value)+((2.*k1*k2)/self.pH.value^2))+(kb*self.boron)/(kb+self.pH.value)+(kw/self.pH.value)-self.pH.value;
                    case "alkalinity"
                        self.co2 = (self.alkalinity-((kb*self.boron)/(kb+self.pH.value))-(kw/self.pH.value)+self.pH.value)/((k1/self.pH.value)+2*((k1*k2)/self.pH.value.^2));
                        self.dic = self.co2*(1+(k1/self.pH.value)+k1*(k2/self.pH.value^2));
                        self.hco3 = self.dic/(1+(self.pH.value/k1)+(k2/self.pH.value));
                        self.co3 = self.dic/(1+(self.pH.value/k2)+((self.pH.value^2)/(k1*k2)));
                    case "dic"
                        self.co2 = self.dic/(1+(k1/self.pH.value)+((k1*k2)/self.pH.value.^2));
                        self.hco3 = self.dic/(1+(self.pH.value/k1)+(k2/self.pH.value));
                        self.co3 = self.dic/(1+(self.pH.value/k2)+((self.pH.value^2)/(k1*k2)));
                        self.alkalinity = self.co2*((k1/self.pH.value)+((2.*k1*k2)/self.pH.value^2))+(kb*self.boron)/(kb+self.pH.value)+(kw/self.pH.value)-self.pH.value;
                end
                self.calculate_CO2();
                self.saturation_state = (self.calcium*self.co3)/self.equilibrium_coefficients.kc.value;
            else
                error("Not implemented yet");
            end
        end
    end
end