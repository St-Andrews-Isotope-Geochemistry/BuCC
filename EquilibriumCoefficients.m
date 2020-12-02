classdef EquilibriumCoefficients < handle
    properties
        k0 = EquilibriumCoefficient();
        k1 = EquilibriumCoefficient();
        k2 = EquilibriumCoefficient();
        kb = EquilibriumCoefficient();
        kw = EquilibriumCoefficient();
        kc = EquilibriumCoefficient();
        ka = EquilibriumCoefficient();
        ks = EquilibriumCoefficient();
        
        kp1 = EquilibriumCoefficient();
        kp2 = EquilibriumCoefficient();
        kp3 = EquilibriumCoefficient();
        
        conditions = Conditions();
    end
    properties (Dependent=true)
        temperature
        salinity
        pressure
        calcium
        magnesium
    end
    properties (Hidden=true)
        property_names = ["k0","k1","k2","kb","kw","kc","ka","ks","kp1","kp2","kp3"];
        calculated = false;
    end
    methods
        % Constructor
        function self = EquilibriumCoefficients()
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
        function set.conditions(self,value);
            self.conditions = value;
            self.setAll("conditions",self.conditions);
        end
        
        
        function set_pressure_correction(self)
            try
                co2_systematics_search = what("CO2_Systematics");
                current_directory = pwd;
                co2_systematics_search = strrep(strrep(co2_systematics_search.path,current_directory,"."),"\","/");
            
                
                raw_file_contents = fileread(co2_systematics_search+"/Configuration/equilibrium_coefficient_pressure_correction.json");                
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
        
        function calculate(self)
            MyAMI_search = what("MyAMI");
            current_directory = pwd;
            if ~isempty(MyAMI_search)
                MyAMI_relative = strrep(strrep(MyAMI_search.path,current_directory,"."),"\","/");
            
                mgca_unit_normalisation = 10^self.conditions.mgca_units_value;
                [k_values,~] = self.run_MyAMI(MyAMI_relative,self.temperature,self.salinity,self.calcium/mgca_unit_normalisation,self.magnesium/mgca_unit_normalisation);
                
                self.kw.value = k_values(1);
                self.k1.value = k_values(2);
                self.k2.value = k_values(3);
                self.kc.value = k_values(4);
                self.kb.value = k_values(5);
                self.ka.value = k_values(6);
                self.k0.value = k_values(7);
                self.ks.value = k_values(8);
                
                self.kw.doPressureCorrection();
                self.k1.doPressureCorrection();
                self.k2.doPressureCorrection();
                self.kc.doPressureCorrection();
                self.kb.doPressureCorrection();
                self.ka.doPressureCorrection();
                self.k0.doPressureCorrection();
                self.ks.doPressureCorrection();
                
                self.calculated = true;
            else
                error("Can't find MyAMI");
            end
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