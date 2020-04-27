within CrudeDistillation.Gamma;

function PartialGamma

input Real a,x ;
output Real ans;
 
algorithm

ans := Modelica.Math.Nonlinear.quadratureLobatto(function Gamma.InputFun(a=a), x, 1000) ;

end PartialGamma;
