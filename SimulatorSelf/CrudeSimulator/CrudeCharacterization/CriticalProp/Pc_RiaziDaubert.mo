within CrudeSimulator.CrudeCharacterization.CriticalProp;
 
function Pc_RiaziDaubert

  input Real PEMe; // Median boiling point of the fraction in K
  input Real d15; //Specific gravity of the fraction
  
  output Real Pc; //Critical Pressure of the fraction
  
  protected
  Real t1;
  
  algorithm
  
  t1 := -0.008505 * PEMe - 4.8014 * d15 + 0.005749 * PEMe * d15;
  Pc := 31958000000.0 * exp(t1) * PEMe ^( -0.4844 )* d15 ^ 4.0846;
            
end Pc_RiaziDaubert;
