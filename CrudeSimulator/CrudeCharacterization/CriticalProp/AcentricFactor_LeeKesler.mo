within CrudeSimulator.CrudeCharacterization.CriticalProp;

function AcentricFactor_LeeKesler

  input Real Tc; //Critical Temerature of the fraction 
  input Real Pc; //Critical Pressure of the fraction
  input Real PEMM; // Molar Mean boiling point of the fraction

  output Real AF; // Accentric factor of the fraction

  algorithm

    AF := (-log(Pc / 101325) - 5.92714 + 6.09648 / (PEMM / Tc) + 1.28862 * log(PEMM / Tc) - 0.169347 * (PEMM / Tc) ^ 6) / (15.2518 - 15.6875 / (PEMM / Tc) - 13.4721 * log(PEMM / Tc) + 0.43577 * (PEMM / Tc) ^ 6);

end AcentricFactor_LeeKesler;
