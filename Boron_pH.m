classdef Boron_pH<handle&Collater&matlab.mixin.Copyable
    % Boron_pH calculates the unknown parameter from d11B_4, d11B_sw, pKb, pH and alpha
    %
    % Boron_pH Properties:
    %   d11B_4 - Delta 11 Borate
    %   d11B_sw - Delta 11 Boron of seawater
    %   epsilon - Fractionation between boric acid and borate, defaults to 27.2‰
    %   pKb - Equivalence point of boric acid and borate, defaults to 8.6
    %   pH - Concentration of hydrogen ions
    % Boron_pH Properties (Dependent):
    %   H - Hydrogen ion concentration (provides access to pH object)
    %   Kb - Provides direct access to pKb.value
    %   alpha - epsilon as a ratio
    %   d11B_3 - Delta 11 Boric acid
    %
    % Boron_pH Methods:
    %   calculate - Determines the missing parameter and calculates its value. Produces an error if there are not enough known parameters, and a warning if all parameters are known. Defaults to using the fully accurate method, specify input of 1 to use simplified method.
    %   
    % Boron_pH Origin:
    %   Written by - Ross Whiteford 1st October 2020
    %   Affiliation - University of St Andrews
    %   Contact - rdmw1@st-andrews.ac.uk
    %   Licensing - Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0), https://creativecommons.org/licenses/by-nc-sa/4.0/
    properties        
        d11B_4 = delta.empty();
        d11B_sw = delta.empty();
        
        epsilon = 27.2;
        
        pKb = pX.empty();
        pH = pX.empty();
         
%         temperature = NaN;
%         salinity = NaN;
%         pressure = NaN;
%         Ca = NaN;
%         Mg = NaN;
    end
    properties (Dependent)
        H;
        Kb;
        alpha;
        d11B_3;
    end
    properties (Hidden)
        validated = 0
    end
    methods
        % Constructor
        function self = Boron_pH(self)
            self.d11B_4 = delta("Boron",NaN);
            self.d11B_sw = delta("Boron",39.61);
            self.pKb = pX(8.6);
            self.pH = pX(NaN);
        end
        
        % Setters
        function set.H(self,value)
            self.pH = -log10(value);
        end
        function set.pH(self,value)
            if isa(value,"pX")
                self.pH = value;
            else
                self.pH = pX(value);
            end
        end
        function set.alpha(self,value)
            self.epsilon = (value-1)*1000;
        end
        function set.d11B_3(self,value)
            self.d11B_4 = delta("B",value-self.epsilon);
        end
        
        % Getters
        function value = get.H(self)
            value = self.pH.value;
        end
        function value = get.Kb(self)
            value = self.pKb.value;
        end
        function value = get.alpha(self)
             value = 1+(self.epsilon/1000);  
        end
        function value = get.d11B_3(self)
            value = delta("B",self.d11B_4.value+self.epsilon);
        end
        
        % d11B_4
        function self = calculate_d11B4_FifthElement(self)
            % Taken from Chapter 5, Boron: The Fifth Element, supplementary information (also on page 109)
            % Should be the same as the full solution (from calculateB4_11_10)
            self.d11B_4.ratio = ((self.H.^2*self.d11B_sw.ratio^2 + 2*self.H.^2*self.d11B_sw.ratio*self.alpha + self.H.^2*self.alpha^2 + 2*self.H.*self.Kb*self.d11B_sw.ratio^2*self.alpha - 2*self.H.*self.Kb*self.d11B_sw.ratio*self.alpha^2 + 8*self.H.*self.Kb*self.d11B_sw.ratio*self.alpha - 2*self.H.*self.Kb*self.d11B_sw.ratio + 2*self.H.*self.Kb*self.alpha + self.Kb.^2*self.d11B_sw.ratio^2*self.alpha^2 + 2*self.Kb.^2*self.d11B_sw.ratio*self.alpha + self.Kb.^2).^(1/2) - self.H.*self.alpha - self.Kb + self.H.*self.d11B_sw.ratio + self.Kb*self.d11B_sw.ratio*self.alpha)./(2*self.alpha*(self.H + self.Kb));           
        end
        function self = calculate_d11B4_simple(self)
            self.d11B_4.ratio = (self.d11B_sw.ratio)./(((self.alpha*self.H)+self.pKb.value)./(self.H+self.pKb.value));
            self.d11B_3.value = self.d11B_4.value+self.epsilon;
        end
        function self = calculate_d11B4(self)
            top1 = 4.*self.alpha.*self.d11B_sw.ratio.*((self.H+self.Kb).^2);
            top2 = (self.H.*(self.d11B_sw.ratio-self.alpha) + (self.Kb.*(self.alpha.*self.d11B_sw.ratio - 1))).^2;
            top3 = (self.H.*(self.d11B_sw.ratio-self.alpha))+(self.Kb*((self.alpha.*self.d11B_sw.ratio)-1));
            bottom = (2.*self.alpha.*(self.H+self.Kb));
            
            B4_11_10 = [((top3+sqrt(top1+top2))./bottom),...
                        ((top3-sqrt(top1+top2))./bottom)];                    
            
            self.d11B_4.ratio = max(B4_11_10);            
        end
        
        % d11B_sw
        function self = calculate_d11B_sw_simple(self)
            self.d11B_sw.ratio = (self.d11B_4.ratio).*(((self.alpha.*self.H)+self.Kb)./(self.H+self.Kb));
        end
        function self = calculate_d11B_sw(self)
            self.d11B_sw.ratio = (((self.alpha.*self.d11B_4.ratio.^2).*(self.Kb+self.H))+(self.d11B_4.ratio.*(self.Kb+(self.alpha.*self.H))))./((self.H.*(self.d11B_4.ratio+1))+(self.Kb.*(self.alpha.*self.d11B_4.ratio+1)));
        end
        
        % pH
        function self = calculate_pH_simple(self)
            self.pH = self.pKb.pValue-log10(-((self.d11B_4.ratio-self.d11B_sw.ratio)./((self.alpha.*self.d11B_4.ratio)-self.d11B_sw.ratio)));
        end
        function self = calculate_pH(self)
            self.pH = self.pKb.pValue-log10( (((self.alpha.*self.d11B_4.ratio)+1).*(self.d11B_4.ratio-self.d11B_sw.ratio))./((self.d11B_4.ratio+1).*(self.d11B_sw.ratio-(self.alpha.*self.d11B_4.ratio))) );
        end
        
        % pKb
        function self = calculate_pKb_simple(self)
            self.pKb = self.pH+log10(-(self.d11B_4.ratio-self.d11B_sw.ratio)./((self.alpha.*self.d11B_4.ratio)-self.d11B_sw.ratio));
        end
        function self = calculate_pKb(self)
            self.pKb = self.pH+log10(((self.alpha.*self.d11B_4.ratio+1).*(self.d11B_4.ratio-self.d11B_sw.ratio))./((self.d11B_4.ratio+1).*(self.d11B_sw.ratio-(self.alpha.*self.d11B_4.ratio))));
        end
        
        % Alpha
        function self = calculate_alpha_simple(self)
            self.alpha = (((self.d11B_sw.ratio.*(self.H+self.Kb))./self.d11B_4.ratio)-self.Kb)./self.H;
        end        
        function self = calculate_alpha(self)
%             one = self.H.*self.d11B_sw.ratio.*(self.d11B_4.ratio+1);
%             two = self.Kb.*(self.d11B_sw.ratio-self.d11B_4.ratio);
%             three = (self.d11B_4.ratio.^2).*(self.H+self.Kb);
%             four = (self.d11B_4.ratio).*(self.H-(self.d11B_sw.ratio.*self.Kb));
            
            self.alpha = ((self.H.*self.d11B_sw.ratio.*(self.d11B_4.ratio+1))+(self.Kb.*(self.d11B_sw.ratio-self.d11B_4.ratio)))./(((self.d11B_4.ratio.^2).*(self.H+self.Kb))+((self.d11B_4.ratio).*(self.H-(self.d11B_sw.ratio.*self.Kb))));
        end
        
        % Validation
        function check_known_values(self)
            number_unknown = double(isnan(self.d11B_4.value))+double(isnan(self.d11B_sw.value))+double(isnan(self.pH.pValue))+double(isnan(self.pKb.pValue))+double(isnan(self.alpha));
            number_known = 5-number_unknown;
            if number_unknown>1
                self.validated = 0;
                error("System underdetermined");
            elseif number_unknown<1
                warning("System overdetermined - skipping");
                self.validated = 2;
            else
                self.validated = 1;
            end
        end
        
        % Combined
        function self = calculate(self,simple)
            %   calculate - Determines the missing parameter and calculates its value. Produces an error if there are not enough known parameters, and a warning if all parameters are known. Defaults to using the fully accurate method, specify input of 1 to use simplified method.
            if nargin<2
                simple = 0;
            end
            for index = 1:numel(self)
                self(index).check_known_values();
                if self(index).validated==1
                    if isnan(self(index).d11B_4.value)
                        if simple
                            self(index).calculate_d11B4_simple();
                        else
                            self(index).calculate_d11B4();
                        end
                    elseif isnan(self(index).d11B_sw.value)
                        if simple
                            self(index).calculate_d11B_sw_simple();
                        else
                            self(index).calculate_d11B_sw();
                        end
                    elseif isnan(self(index).pH.pValue)
                        if simple
                            self(index).calculate_pH_simple();
                        else
                            self(index).calculate_pH();
                        end
                    elseif isnan(self(index).pKb.pValue)
                        if simple
                            self(index).calculate_pKb_simple();
                        else
                            self(index).calculate_pKb();
                        end
                    elseif isnan(self(index).alpha)
                        if simple
                            self(index).calculate_alpha_simple();
                        else
                            self(index).calculate_alpha();
                        end
                    else
                        error("Unforeseen error...");
                    end
                end
            end
        end
    end
end