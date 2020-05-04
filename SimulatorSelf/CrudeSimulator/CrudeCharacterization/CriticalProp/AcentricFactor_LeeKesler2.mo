within CrudeSimulator.CrudeCharacterization.CriticalProp;

function AcentricFactor_LeeKesler2 // http://dwsim.inforside.com.br/docs/mobile/help/thermo.htm

input Real Tc;
input Real Pc;
input Real Tb;

output Real AF;

protected

Real num, de, Tbr;

algorithm

Tbr := Tb/Tc ;

num := -log(Pc/101325) - 5.92714 + 6.0948/Tbr + 1.28862*log(Tbr) - 0.169347*Tbr^6 ;
de := 15.2518 - 15.6875/Tbr - 13.472*log(Tbr) + 0.43577*Tbr^6 ;

AF := num/de ;

end AcentricFactor_LeeKesler2;
