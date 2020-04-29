within Pseudocomponents;

model GeneralProperties

parameter Integer SN "Serial Number" ;
parameter Real Tc (unit = "K" ) "Critical Temperature" ;
parameter Real Pc (unit = "Pa" ) "Critical Pressure" ;
parameter Real Tb (unit = "K" ) "Boiling Temperature" ;
parameter Real MW (unit = "kg/kmol") "Molecular Weight" ;
parameter Real AF (unit = "-") "Acentric Factor" ;
parameter Real VP[6] (unit = "Pa" ) "Vapour Pressure constants" ;
parameter Real HOV[6] (unit = "J/kmol") "Heat of vaparization constants" ;
parameter Real VapCp[6] (unit = "J/kmol.K") "Heat capacity of Vapour constants " ;
parameter Real LiqCp[6] (unit = "J/kmol.K") "Heat capacity of Liquid constants " ;

end GeneralProperties;
