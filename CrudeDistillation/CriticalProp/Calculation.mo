within CrudeDistillation.CriticalProp;

model Calculation

extends BulkCalculation.Raizi;
extends input_values;

Real Pc[n], d15[n], Tc[n], AF[n] ;

equation

for i in 1:n loop
  
  d15[i] = d15_Riazi( MWai[i] ) ;
  Tc[i] = Tc_Riazi( Tbai[i], d15[i] ) ;
  Pc[i] = Pc_RiaziDaubert( Tbai[i], d15[i] ) ;
  AF[i] = AcentricFactor_LeeKesler( Tc[i], Pc[i], Tbai[i] ) ;
  
end for;
  
end Calculation;
