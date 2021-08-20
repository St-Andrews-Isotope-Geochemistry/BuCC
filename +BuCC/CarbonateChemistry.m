classdef CarbonateChemistry < handle&Geochemistry_Helpers.Collator
    properties
        pH
        
        dic = NaN;
        alkalinity = NaN;
        
        oceanic_co2 = NaN;
        hco3 = NaN;
        co3 = NaN;
        atmospheric_co2 = NaN;
        
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
        oceanic_pressure
        atmospheric_pressure
        calcium
        magnesium
        calculable
    end
    methods
        % Constructor
        function self = CarbonateChemistry(varargin)
            parser = inputParser;
            properties = ["temperature","salinity","oceanic_pressure","atmospheric_pressure","calcium","magnesium","conditions","dic","alkalinity","pH","oceanic_co2","hco3","co3","atmospheric_co2","saturation_state"];
            
            for property = properties
                addOptional(parser,property,NaN);
            end
            
            parse(parser,varargin{:});
            
            for property = properties
                self.(property) = parser.Results.(property);
            end
            
            if isnan(parser.Results.conditions)
                self.conditions = BuCC.Conditions(varargin{:});
            end
            
            self.equilibrium_coefficients = BuCC.EquilibriumCoefficients(true);
            json_pressure = jsondecode(fileread("equilibrium_coefficient_pressure_correction.json"));
            json_function = jsondecode(fileread("equilibrium_coefficient_functions.json"));
            for property_name = self.equilibrium_coefficients.property_names
                self.equilibrium_coefficients.(property_name).parsePressureCorrectionsAndFunctions(json_pressure,json_function);
            end

            self.equilibrium_coefficients.valid_file_found = true;
            self.equilibrium_coefficients.conditions = self.conditions;
            
            self.atmospheric_co2 = BuCC.AtmosphericCO2();
            self.atmospheric_co2.conditions = self.conditions;
            
            self.pH = Geochemistry_Helpers.pX();
        end
        
        % Setters
        function set.temperature(self,value)
            self.conditions.temperature = value;
        end
        function set.salinity(self,value)
            self.conditions.salinity = value;
        end
        function set.oceanic_pressure(self,value)
            self.conditions.oceanic_pressure = value;
        end
        function set.atmospheric_pressure(self,value)
            self.conditions.atmospheric_pressure = value;
        end
        function set.calcium(self,value)
            self.conditions.calcium = value;
        end
        function set.magnesium(self,value)
            self.conditions.magnesium = value;
        end
        function set.units(self,value)
            self.units = value;
            self.getUnitsValue();
            self.units_set = true;
        end
        function set.pH(self,value)
            if isnumeric(value)
                self.pH.pValue = value;
            else
                self.pH = value;
            end
        end
        
        function setConditions(self,conditions)
            self.conditions.temperature = conditions.temperature;
            self.conditions.salinity = conditions.salinity;
            self.conditions.oceanic_pressure = conditions.oceanic_pressure;
            self.conditions.calcium = conditions.calcium;
            self.conditions.magnesium = conditions.magnesium;
            
            self.conditions.atmospheric_pressure = conditions.atmospheric_pressure;
        end
        
         % Getters
        function output = get.temperature(self)
            output = self.conditions.temperature;
        end
        function output = get.salinity(self)
            output = self.conditions.salinity;
        end
        function output = get.oceanic_pressure(self)
            output = self.conditions.oceanic_pressure;
        end
        function output = get.atmospheric_pressure(self)
            output = self.conditions.atmospheric_pressure;
        end
        function output = get.calcium(self)
            output = self.conditions.calcium;
        end
        function output = get.magnesium(self)
            output = self.conditions.magnesium;
        end
        function output = get.calculable(self)
            known_properties = self.getKnownProperties();
            if numel(known_properties)==2
                output = true;
            else
                output = false;
            end
        end
        
        function known_properties = getKnownProperties(self)
            known_properties = [];
            
            if ~isnan(self.pH.value)
                known_properties = [known_properties,"pH"];
            end
            if ~isnan(self.dic)
                known_properties = [known_properties,"dic"];
            end
            if ~isnan(self.alkalinity)
                known_properties = [known_properties,"alkalinity"];
            end
            if ~isnan(self.oceanic_co2) || ~isnan(self.atmospheric_co2)
                known_properties = [known_properties,"oceanic_co2"];
            end
            if ~isnan(self.hco3)
                known_properties = [known_properties,"hco3"];
            end
            if ~isnan(self.co3)
                known_properties = [known_properties,"co3"];
            end
            if ~isnan(self.saturation_state)
                known_properties = [known_properties,"saturation_state"];
            end            
        end
        
        % Units
        function self = getUnitsValue(self)
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
        function self = estimateUnits(self,parameter_1,parameter_2)
            if ~self.units_set
                units_char = char(self.units);
                typical_values_map = containers.Map(["dic","alkalinity","oceanic_co2","hco3","co3","atmospheric_co2"],[2000e-6,2300e-6,10e-6,1800e-6,200e-6,300e-6]);
                if units_char(1)=="x"
                    
                    typical_order_of_magnitude = 3*round(floor(log10(typical_values_map(lower(parameter_1))))/3,0);
                    if ~strcmp(parameter_1,"atmospheric_co2")
                        input_order_of_magnitude = 3*round(floor(log10(self.(parameter_1)))/3,0);
                    else
                        input_order_of_magnitude = 3*round(floor(log10(self.atmospheric_co2.mole_fraction))/3,0);
                    end
                    difference_order_of_magnitude = input_order_of_magnitude-typical_order_of_magnitude;
                    self.units_value = difference_order_of_magnitude;

                    if difference_order_of_magnitude==6
                        units_char(1) = "μ";
                    elseif difference_order_of_magnitude==3
                        units_char(1) = "m";
                    elseif difference_order_of_magnitude==0
                        units_char(1) = " ";
                    end
                    self.units = units_char;
                end

                if nargin==3                    
                    typical_order_of_magnitude = 3*round(floor(log10(typical_values_map(parameter_2)))/3,0);
                    input_order_of_magnitude = 3*round(floor(log10(self.(parameter_2)))/3,0);
                    difference_order_of_magnitude = input_order_of_magnitude-typical_order_of_magnitude;

                    if difference_order_of_magnitude==6
                        assert(units_char(1)=="μ","Units mismatch");
                    elseif difference_order_of_magnitude==3
                        assert(units_char(1)=="m","Units mismatch");
                    elseif difference_order_of_magnitude==0
                        assert(units_char(1)==" ","Units mismatch");
                    end
                end
            end
        end
        
        % Calculations
        function self = calculateBoronFromSalinity(self)
            self.boron = (0.0002414/10.811)*(self.salinity/1.80655);
        end
        function self = calculate_CO2(self)
            for self_index=1:numel(self)
                if ~isnan(self(self_index).atmospheric_co2)
                    self(self_index).estimateUnits("atmospheric_co2");
                elseif ~isnan(self(self_index).oceanic_co2)
                    self(self_index).estimateUnits("oceanic_co2");
                end
                
                
                k0 = self(self_index).equilibrium_coefficients.k0.value;
                unit_normalisation = 10^self(self_index).units_value;
                
                atmospheric_co2 = self(self_index).atmospheric_co2.partial_pressure/unit_normalisation;
                ocean_co2 = self(self_index).oceanic_co2/unit_normalisation;
                
                if isnan(ocean_co2) && ~isnan(atmospheric_co2)
                    ocean_co2 = atmospheric_co2*k0;
                elseif ~isnan(self(self_index).oceanic_co2) && isnan(self(self_index).atmospheric_co2)
                    atmospheric_co2 = (ocean_co2/k0);
                elseif ~isnan(self(self_index).oceanic_co2) && ~isnan(self(self_index).atmospheric_co2)
                    atmospheric_co2_guess = ocean_co2/k0;
                    if abs(atmospheric_co2_guess-atmospheric_co2)>1e-12
                        error("Inconsistent ocean CO2 and atmospheric CO2");
                    end
                end
                self(self_index).oceanic_co2 = ocean_co2*unit_normalisation;
                self(self_index).atmospheric_co2.partial_pressure = atmospheric_co2*unit_normalisation;
            end
        end
        function self = calculate(self)
            for self_index=1:numel(self)
                mgca_unit_normalisation = 10^self(self_index).conditions.mgca_units_value;
                calcium = self(self_index).calcium/mgca_unit_normalisation;
                magnesium = self(self_index).magnesium/mgca_unit_normalisation;
                
                if ~self(self_index).equilibrium_coefficients.calculated
                    self(self_index).equilibrium_coefficients.calculate();
                end
                known_properties = self(self_index).getKnownProperties();

                if isnan(self(self_index).boron)
                    self(self_index).calculateBoronFromSalinity();
                end
                if ~isnan(self(self_index).atmospheric_co2)
                    self(self_index).atmospheric_co2.calculate();
                end
                if ~isnan(self(self_index).oceanic_co2) || ~isnan(self(self_index).atmospheric_co2)
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
                        case "oceanic_co2"
                            self(self_index).estimateUnits("oceanic_co2");
                            unit_normalisation = 10^self(self_index).units_value;

                            oceanic_co2 = (self(self_index).oceanic_co2/unit_normalisation);

                            dic = oceanic_co2*(1+(k1/pH)+k1*(k2/pH^2));
                            hco3 = dic/(1+(pH/k1)+(k2/pH));
                            co3 = dic/(1+(pH/k2)+((pH^2)/(k1*k2)));
                            alkalinity = oceanic_co2*((k1/pH)+((2*k1*k2)/pH^2))+(kb*self(self_index).boron)/(kb+pH)+(kw/pH)-pH;
                        case "hco3"
                            self(self_index).estimateUnits("hco3");
                            unit_normalisation = 10^self(self_index).units_value;

                            hco3 = (self(self_index).hco3/unit_normalisation);

                            dic = hco3*(1+(pH/k1)+(k2/pH));
                            oceanic_co2 = dic/(1+(k1/pH)+((k1*k2)/pH^2));
                            co3 = dic/(1+(pH/k2)+((pH^2)/(k1*k2)));
                            alkalinity = oceanic_co2*((k1/pH)+((2*k1*k2)/pH^2))+(kb*self(self_index).boron)/(kb+pH)+(kw/pH)-pH;
                        case "co3"
                            self(self_index).estimateUnits("co3");
                            unit_normalisation = 10^self(self_index).units_value;

                            co3 = (self(self_index).co3/unit_normalisation);

                            dic = co3*(1+(pH/k2)+(pH^2/(k1*k2)));
                            hco3 = dic/(1+(pH/k1)+(k2/pH));
                            oceanic_co2 = dic/(1+(k1/pH)+((k1*k2)/pH^2));
                            alkalinity = oceanic_co2*((k1/pH)+((2*k1*k2)/pH^2))+(kb*self(self_index).boron)/(kb+pH)+(kw/pH)-pH;
                        case "alkalinity"
                            self(self_index).estimateUnits("alkalinity");
                            unit_normalisation = 10^self(self_index).units_value;

                            alkalinity = (self(self_index).alkalinity/unit_normalisation);

                            oceanic_co2 = (alkalinity-((kb*self(self_index).boron)/(kb+pH))-(kw/pH)+pH)/((k1/pH)+2*((k1*k2)/pH^2));
                            dic = oceanic_co2*(1+(k1/pH)+k1*(k2/pH^2));
                            hco3 = dic/(1+(pH/k1)+(k2/pH));
                            co3 = dic/(1+(pH/k2)+((pH^2)/(k1*k2)));
                        case "dic"
                            self(self_index).estimateUnits("dic");
                            unit_normalisation = 10^self(self_index).units_value;

                            dic = (self(self_index).dic/unit_normalisation);

                            oceanic_co2 = dic/(1+(k1/pH)+((k1*k2)/pH^2));
                            hco3 = dic/(1+(pH/k1)+(k2/pH));
                            co3 = dic/(1+(pH/k2)+((pH^2)/(k1*k2)));
                            alkalinity = oceanic_co2*((k1/pH)+((2*k1*k2)/pH^2))+((kb*self(self_index).boron)/(kb+pH))+(kw/pH)-pH;
                        otherwise
                            error("Not implemented yet");
                    end
                    self(self_index).dic = dic*unit_normalisation;
                    self(self_index).alkalinity = alkalinity*unit_normalisation;
                    self(self_index).oceanic_co2 = oceanic_co2*unit_normalisation;
                    self(self_index).hco3 = hco3*unit_normalisation;
                    self(self_index).co3 = co3*unit_normalisation;

                    self(self_index).calculate_CO2();
                    self(self_index).saturation_state = (calcium*co3)/self(self_index).equilibrium_coefficients.kc.value;
                elseif ~isnan(self(self_index).dic) && ~isnan(self(self_index).oceanic_co2)
                    self(self_index).estimateUnits("dic");
                    unit_normalisation = 10^self(self_index).units_value;
                            
                    dic = (self(self_index).dic/unit_normalisation);
                    oceanic_co2 = (self(self_index).oceanic_co2/unit_normalisation);                    
                    
                    p2 = dic-oceanic_co2;
                    p1 = -oceanic_co2*k1;
                    p0 = -oceanic_co2*k1*k2;
                    p = [p2,p1,p0];
                    r = roots(p);
                    pH = max(real(r));
                    
                    dic = oceanic_co2*(1+(k1/pH)+k1*(k2/pH^2));
                    hco3 = dic/(1+(pH/k1)+(k2/pH));
                    co3 = dic/(1+(pH/k2)+((pH^2)/(k1*k2)));
                    alkalinity = oceanic_co2*((k1/pH)+((2*k1*k2)/pH^2))+(kb*self(self_index).boron)/(kb+pH)+(kw/pH)-pH;
                        
                    self(self_index).dic = dic*unit_normalisation;
                    self(self_index).alkalinity = alkalinity*unit_normalisation;
                    self(self_index).oceanic_co2 = oceanic_co2*unit_normalisation;
                    self(self_index).hco3 = hco3*unit_normalisation;
                    self(self_index).co3 = co3*unit_normalisation;

                    self(self_index).saturation_state = ((calcium*co3)/self(self_index).equilibrium_coefficients.kc.value);
                    self(self_index).pH.value = pH;
                elseif ~isnan(self(self_index).alkalinity) && ~isnan(self(self_index).oceanic_co2)
                    self(self_index).estimateUnits("alkalinity");
                    unit_normalisation = 10^self(self_index).units_value;
                    
                    alkalinity = (self(self_index).alkalinity/unit_normalisation);
                    oceanic_co2 = (self(self_index).oceanic_co2/unit_normalisation);
                    
                    p4 = 1.;
                    p3 = kb+alkalinity;
                    p2 = alkalinity*kb-oceanic_co2*k1-kb*self.boron-kw;
                    p1 = -oceanic_co2*kb*k1-oceanic_co2*2.*k1*k2-kw*kb;
                    p0 = -2.*oceanic_co2*kb*k1*k2;
                    p = [p4,p3,p2,p1,p0];
                    r = roots(p);
                    pH = max(real(r));
                    
                    dic = oceanic_co2*(1+(k1/pH)+k1*(k2/pH^2));
                    hco3 = dic/(1+(pH/k1)+(k2/pH));
                    co3 = dic/(1+(pH/k2)+((pH^2)/(k1*k2)));
                    
                    self(self_index).dic = dic*unit_normalisation;
                    self(self_index).alkalinity = alkalinity*unit_normalisation;
                    self(self_index).oceanic_co2 = oceanic_co2*unit_normalisation;
                    self(self_index).hco3 = hco3*unit_normalisation;
                    self(self_index).co3 = co3*unit_normalisation;
                    
                    self(self_index).saturation_state = ((calcium*co3)/self(self_index).equilibrium_coefficients.kc.value);
                    self(self_index).pH.value = pH;
                elseif ~isnan(self(self_index).oceanic_co2) && ~isnan(self(self_index).co3)
                    self(self_index).estimateUnits("co3","oceanic_co2");
                    unit_normalisation = 10^self(self_index).units_value;
                    
                    co3 = (self(self_index).co3/unit_normalisation);
                    oceanic_co2 = (self(self_index).oceanic_co2/unit_normalisation);                    
                    
                    p4 = -co3/k1/k2;
                    p3 = -co3/k2;
                    p2 = oceanic_co2-co3;
                    p1 = oceanic_co2*k1;
                    p0 = oceanic_co2*k1*k2;
                    p = [p4,p3,p2,p1,p0];
                    r = roots(p);
                    pH = max(real(r));
                    
                    dic = oceanic_co2*(1.+k1/pH+k1*k2/pH/pH);
                    hco3 = dic/(1+pH/k1+k2/pH);
                    alkalinity = oceanic_co2*(k1/pH+2.*k1*k2/pH/pH)+kb*self(self_index).boron/(kb+pH)+kw/pH-pH;
                    
                    self(self_index).dic = dic*unit_normalisation;
                    self(self_index).alkalinity = alkalinity*unit_normalisation;
                    self(self_index).oceanic_co2 = oceanic_co2*unit_normalisation;
                    self(self_index).hco3 = hco3*unit_normalisation;
                    self(self_index).co3 = co3*unit_normalisation;
                    self(self_index).saturation_state = ((calcium*co3)/self(self_index).equilibrium_coefficients.kc.value);
                    self(self_index).pH.value = pH;
                else
                    error("Not implemented yet");
                end
                self(self_index).atmospheric_co2.calculate();
            end
        end
        
        % Display
        function show(self,parameter)
            if parameter=="alkalinity" || parameter=="dic" || parameter=="oceanic_co2" || parameter=="hco3" || parameter=="co3"
                disp(join([num2str(self.(parameter)),self.units]));
            elseif parameter=="temperature"
                disp(join([num2str(self.(parameter)),"C"],""));
            else
                disp(self.(parameter));
            end
        end        
    end
end