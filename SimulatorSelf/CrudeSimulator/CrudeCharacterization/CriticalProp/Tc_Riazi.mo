within CrudeSimulator.CrudeCharacterization.CriticalProp;

function Tc_Riazi

  input Real PEMe; // Median boiling point of the fraction in K
  input Real d15; //Specific gravity of the fraction
  
  output Real Tc; //Critical Temperature of the fraction
  
  protected
  Real t1;
  
  algorithm
  
  t1 := -0.00069 * PEMe - 1.4442 * d15 + 0.000491 * PEMe * d15;
  Tc := 35.9413 * exp(t1) * PEMe ^ 0.7293 * d15 ^ 1.2771;
            

end Tc_Riazi;
