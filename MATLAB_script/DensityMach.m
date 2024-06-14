%/////////////////////////////////////////////////////////////////////////%
%                                                                         %
%   - Name : Air Density function                                         %
%   - Compute Air Density at Specific Altitude                            %
%                                                                         %
%   - Inputs :  (1) Alt = Altitude (meter)                                %
%                                                                         %
%   - Outputs : (1) Rho = Air Density (kg/m^3)                            %
%                                                                         %
%                                   2015. 02. 05. Created by Hong, S. M.  %
%                                                                         %
%/////////////////////////////////////////////////////////////////////////%

%.. Environment Function

function [Rho Mach] = DensityMach( Alt, Vel  )

%.. Declare Global Variables

% global          R              P_sl

gamma           =       1.4     ;
R               =       287.058 ;                                          % Ideal Gas Constant of Air  (m^2/(K*s^2))    R_air = R_universal / M_air 
P_sl            =       101325  ;                                          % Pressure at Sea Level      (Pa)

%.. Compute Environment Variables

if  Alt >= 0
    
    if  Alt < 11000                                                        % Troposphere Case
        
        T           =       288.15 - 0.0065 * Alt ;                        % Temperature                (K)
        Pressure    =       P_sl * ( T / 288.15 ) ^ 5.2559 ;               % Atmospheric Pressure       (Pa)
        
    else                                                                   % Stratosphere Case
         
        T           =       216 ;                                          % Temperature                (K)
        Pressure    =       22630 * exp( - 0.00015769 * ( Alt - 11000 ) ) ;% Atmospheric Pressure       (Pa)
        
    end
    
    Rho             =       Pressure / ( R * T )  ;                        % Atmospheric Density        (kg/m^3)
    
else
    
    T               =       288.15  ;
    Pressure        =       P_sl    ;
    Rho             =       Pressure / ( R * T ) ;
    
end

Sound               = sqrt(gamma * R * T);
Mach                = Vel/Sound          ;

