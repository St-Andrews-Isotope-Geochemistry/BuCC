classdef EquilibriumCoefficients < handle&Geochemistry_Helpers.Collator
    properties
        k0
        k1
        k2
        kb
        kw
        kc
        ka
        ks
        kf
        
        kp1
        kp2
        kp3
        
        MyAMI;
        
        sulphate = 0.02824;
        fluoride = 7e-5;
        
        conditions;
    end
    properties (Dependent=true)
        temperature
        salinity
        pressure
        calcium
        magnesium
    end
    properties (Hidden=true)
        property_names = ["k0","k1","k2","kb","kw","kc","ka","ks","kp1","kp2","kp3","kf"];
        calculated = false;
    end
    methods
        % Constructor
        function self = EquilibriumCoefficients()
            self.conditions = BuCC.Conditions();
            
            self.k0 = BuCC.EquilibriumCoefficient();
            self.k1 = BuCC.EquilibriumCoefficient();
            self.k2 = BuCC.EquilibriumCoefficient();
            self.kb = BuCC.EquilibriumCoefficient();
            self.kw = BuCC.EquilibriumCoefficient();
            self.kc = BuCC.EquilibriumCoefficient();
            self.ka = BuCC.EquilibriumCoefficient();
            self.ks = BuCC.EquilibriumCoefficient();
            self.kf = BuCC.EquilibriumCoefficient();
            
            self.kp1 = BuCC.EquilibriumCoefficient();
            self.kp2 = BuCC.EquilibriumCoefficient();
            self.kp3 = BuCC.EquilibriumCoefficient();
            
            self.set_pressure_correction();
            self.setAll("conditions",self.conditions);
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
        function set.conditions(self,value)
            self.conditions = value;
            self.setAll("conditions",self.conditions);
        end        
        
        function set_pressure_correction(self)
            try
                bucc_search = what("+Bucc");
                current_directory = pwd;
                bucc_search = strrep(strrep(join([bucc_search(1).path,"\Configuration"],""),current_directory,"."),"\","/");
            
                
                raw_file_contents = fileread(bucc_search+"/equilibrium_coefficient_pressure_correction.json");
                json_file_contents = jsondecode(raw_file_contents);
                valid_json_file_found = 1;
            catch
                valid_json_file_found = 0;
            end
            
            if valid_json_file_found==1
                for property_index = 1:numel(self.property_names)
                    self.(self.property_names(property_index)).pressure_correction = json_file_contents.(self.property_names(property_index));
                end
            else
                for property_index = 1:numel(self.property_names)
                    self.(self.property_names(property_index)).pressure_correction = [NaN,NaN,NaN,NaN,NaN];
                end
            end
        end
        function setAll(self,parameter,value)
            for property_index = 1:numel(self.property_names)
                self.(self.property_names(property_index)).(parameter) = value;
            end
        end
        
        function [k_values,k_correction] = getEquilibriumCoefficients(self,method)
            if nargin<2
                method = "MyAMI";
            end
            BuCC_search = what("+BuCC");
            current_directory = pwd; 
            MyAMI_relative = strrep(strrep(join([BuCC_search(1).path,"\+MyAMI"],""),current_directory,"."),"\","/");
            mgca_unit_normalisation = 10^self.conditions.mgca_units_value;
            calcium = self.calcium/mgca_unit_normalisation;
            magnesium = self.magnesium/mgca_unit_normalisation;
            if strcmp(method,"MyAMI")
                [k_values,k_correction] = self.run_MyAMI(MyAMI_relative,self.temperature,self.salinity,calcium,magnesium);
            elseif strcmp(method,"Precalculated")
                try
                    calcium_magnesium = readmatrix(MyAMI_relative+"./Precalculated.xls","Range",join(["A:B"],""));
                    raw_file_contents = fileread(BuCC_search(1).path+"/Configuration/equilibrium_coefficient_functions.json");
                    json_file_contents = jsondecode(raw_file_contents);
                    json_fieldnames = fieldnames(json_file_contents);
                    for fieldname_index = 1:numel(json_fieldnames)
                        self.(json_fieldnames{fieldname_index}).function_handle = str2func(json_file_contents.(json_fieldnames{fieldname_index}));
                    end                    
                catch
                    error("File not found");
                end
                header_rows = sum(isnan(calcium_magnesium(:,1)));
                calcium_known = calcium_magnesium(header_rows+1:end,1);
                magnesium_known = calcium_magnesium(header_rows+1:end,2);
                
                calcium_unique = unique(calcium_known);
                magnesium_unique = unique(magnesium_known);
                
                calcium_resolution = calcium_unique(2)-calcium_unique(1);
                calcium_range = max(calcium_unique)-min(calcium_unique);
                
                magnesium_resolution = magnesium_unique(2)-magnesium_unique(1);
                magnesium_range = max(magnesium_unique)-min(magnesium_unique);
                
                if mod(calcium,calcium_resolution)==0 && mod(magnesium,magnesium_resolution)==0 % There's an exact result in the spreadsheet
                    index = header_rows+(1+magnesium/magnesium_resolution)+((calcium/calcium_resolution)*(1+magnesium_range/magnesium_resolution));
                    coefficients = readmatrix(MyAMI_relative+"./Precalculated.xls","Range",join(["D",num2str(index),":BV",num2str(index)],""));
                else
                    calcium_query = [floor(calcium*1000)/1000,ceil(calcium*1000)/1000];
                    magnesium_query = [floor(magnesium*1000)/1000,ceil(magnesium*1000)/1000];
                    
                    for calcium_query_index = 1:numel(calcium_query)
                        for magnesium_query_index = 1:numel(magnesium_query)
                            index = header_rows+(1+magnesium_query(magnesium_query_index)/magnesium_resolution)+((calcium_query(calcium_query_index)/calcium_resolution)*(1+magnesium_range/magnesium_resolution));
                            coefficients(magnesium_query_index,calcium_query_index,:) = readmatrix(MyAMI_relative+"./Precalculated.xls","Range",join(["D",num2str(index),":BV",num2str(index)],""));
                        end
                    end
                    
                    for coefficient_index = 1:size(coefficients,3)
                        interpolated_coefficient(coefficient_index) = 1/((calcium_query(2)-calcium_query(1))*(magnesium_query(2)-magnesium_query(1))) * [calcium_query(2)-calcium,calcium-calcium_query(1)] * coefficients(:,:,coefficient_index) * [magnesium_query(2)-magnesium;magnesium-magnesium_query(1)];
                    end
                    coefficients = interpolated_coefficient;    
                    
                end
                coefficients = [coefficients,NaN];
                k_order = ["k0","k1","k2","kb","kw","kc","ka","ks"];
                k_count = 1;
                coefficient_start = 1;
                for coefficient_index = 1:numel(coefficients)
                    if ~isnan(coefficients(coefficient_index))
                        continue
                    else
                        self.(k_order(k_count)).function_coefficients = coefficients(coefficient_start:coefficient_index-1);
                        k_count = k_count+1;
                        coefficient_start = coefficient_index+1;
                    end
                end
                
                for ck = k_order
                    self.(ck).calculate();
                end
            end
        end
        
        function calculate(self)
            if ~self.calculated
%                 if ~isempty(BuCC_search)
%                     MyAMI_relative = strrep(strrep(join([BuCC_search(1).path,"\+MyAMI"],""),current_directory,"."),"\","/");
% 
                    mgca_unit_normalisation = 10^self.conditions.mgca_units_value;
%                     [k_values,~] = self.run_MyAMI(MyAMI_relative,self.temperature,self.salinity,self.calcium/mgca_unit_normalisation,self.magnesium/mgca_unit_normalisation);

                    if isempty(self.MyAMI)
                        self.MyAMI = MyAMI.MyAMI();
                    end                        
                    self.MyAMI.calculate(self.conditions.temperature,self.conditions.salinity,self.calcium/mgca_unit_normalisation,self.magnesium/mgca_unit_normalisation,"Precalculated",true);
                    k_values = self.MyAMI.results;
                    
                    self.k0.value = k_values("k0");
                    self.k1.value = k_values("k1");
                    self.k2.value = k_values("k2");
                    self.kw.value = k_values("kw");
                    self.kb.value = k_values("kb");
                    self.kc.value = k_values("kc");
                    self.ka.value = k_values("ka");
                    self.ks.value = k_values("ks");

                    self.kf.value = 0.001764409566690456265466990793;
                    self.kf.doPressureCorrection();

                    self.ks.doPressureCorrection();

                    tb = 1+self.ks.correction+self.ks.value/self.sulphate+(self.sulphate*self.ks.correction)/self.ks.value;
                    t = (self.fluoride/self.kf.value)*((self.ks.value/self.sulphate)*self.kf.correction + self.kf.correction);
                    b = (self.fluoride/self.kf.value)*((self.ks.value/self.sulphate) + self.ks.correction);

                    scale_correction = (tb+t)/(tb+b);

                    self.kw.doScaleCorrectedPressureCorrection(scale_correction);
                    self.k1.doScaleCorrectedPressureCorrection(scale_correction);
                    self.k2.doScaleCorrectedPressureCorrection(scale_correction);
                    self.kc.doScaleCorrectedPressureCorrection(scale_correction);
                    self.kb.doScaleCorrectedPressureCorrection(scale_correction);
                    self.ka.doScaleCorrectedPressureCorrection(scale_correction);
                    self.k0.doScaleCorrectedPressureCorrection(scale_correction);

                    self.calculated = true;
                else
                    error("Can't find MyAMI");
                end
%             end
        end
    end
    methods (Static)
        function [k_values,k_values_correction] = run_MyAMI(MyAMI_path,temperature,salinity,calcium,magnesium);
            command = join(["python ",MyAMI_path,"/PITZER.py ",temperature," ",salinity," ",calcium," ",magnesium],"");
            [status,result] = system(command);
            if status==0
            % Should check status and result for reasonable values
                result_cleaned = erase(string(result(2:end-2)),newline);
                result_values = str2num(result_cleaned);
                input_values = result_values(1:4);
                output_values = result_values(5:end-1); % Last value unused in the examples...
                output_matrix = reshape(output_values,[],3);
                
                k_values_correction = output_matrix(:,1)./output_matrix(:,3);
                k_values = output_matrix(:,2).*k_values_correction;
            else
                error(join(["Problem running Python:",result]," "));
            end
        end
    end
end