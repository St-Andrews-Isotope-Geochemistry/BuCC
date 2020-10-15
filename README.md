# Boron Systematics
Boron Systematics is a Matlab class to take any four of the following:  
- pH (or [H]<sup>+</sup>])  
- pKb (or Kb)  
- &alpha; (or &epsilon;)  
- &delta;<sup>11</sup>B<sub>sw</sub>  
- &delta;<sup>11</sup>B<sub>borate</sub>  

and calculates the remaining parameter.

## Dependencies
The Matlab class ```Boron_pH``` depends on the ```Geochemistry_Helpers``` [repo](https://github.com/St-Andrews-Isotope-Geochemistry/Geochemistry_Helpers).

```Geochemistry_Helpers``` is a series of Matlab classes to make processing Geochemical data easier.

## Getting Started  
The simplest usage of ```Boron_pH``` is as follows:  
Ensure both ```Boron_pH``` and the ```Geochemistry_Helpers``` subdirectory are on the Matlab path.
Then run:
```MATLAB
data = Boron_pH();      % Instantiate an object
data.d11B_4.value = 20; % Fill in either a pH or a d11B_4
data.calculate();       % Call the calculate method
disp(data.pH.pValue)    % Show the result
```

Similarly, the process can be reversed:
```MATLAB
data = Boron_pH();      % Instantiate an object
data.pH.pValue = 8.2;   % Fill in either a pH or a d11B_4
data.calculate();       % Call the calculate method
disp(data.d11B_4.value) % Show the result
```

Alternatively, any of the assumed values can be overwritten e.g.:
```MATLAB
data = Boron_pH();       % Instantiate an object
data.pH.pValue = 8.2;    % Fill in either a pH or a d11B_4
data.d11B_sw.value = 35; % Overwrite assumed d11B_sw
data.calculate();        % Call the calculate method
disp(data.d11B_4.value)  % Show the result
```

Or if one of the assumed values is unknown, overwrite with NaN e.g.:
```MATLAB
data = Boron_pH();        % Instantiate an object
data.pH.pValue = 8.2;     %
data.d11B_4.Value = 20;   % Fill in both pH and d11B_4
data.d11B_sw.value = NaN; % Overwrite assumed d11B_sw
data.calculate();         % Call the calculate method
disp(data.d11B_sw.value)   % Show the result
```

## Multiple Values
If you're dealing with arrays of data rather than single values, then you can make arrays of ```Boron_pH``` objects.
```MATLAB
data = [Boron_pH(),Boron_pH(),Boron_pH()];  % Instantiate an  array of objects

data(1).pH.pValue = 8.0;
data(2).pH.pValue = 8.2;
data(3).pH.pValue = 8.4;                    % Fill in either a pH or a d11B_4

data.calculate();                           % Call the calculate method

disp(data(1).d11B_4.value)
disp(data(2).d11B_4.value)
disp(data(3).d11B_4.value)                  % Show the results
```


## Useful Information
```Boron_pH``` assumes modern day values for pKb (8.6), &epsilon; (27.2), and &delta;<sup>11</sup>B<sub>sw</sub> (39.61). These can be overwritten if required.
