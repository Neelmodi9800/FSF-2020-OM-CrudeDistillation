within CrudeSimulator.CrudeCharacterization.Gamma;

function InputFun
extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction;

input Real a;

algorithm

y := u^(a-1) * exp(-u) ;

end InputFun;
