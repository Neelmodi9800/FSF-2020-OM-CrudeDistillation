within CrudeDistillation.CriticalProp;

function d15_Riazi

input Real MW; // Molecular Weight of Fraction

output Real d15; // Specific gravity at 15 °C

algorithm
 
d15 := 1.07 - exp(3.56073 - 2.93886 * MW ^ 0.1);

end d15_Riazi;
