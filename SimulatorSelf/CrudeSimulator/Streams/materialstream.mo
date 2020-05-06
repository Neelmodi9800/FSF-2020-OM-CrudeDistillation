within CrudeSimulator.Streams;


model materialstream "Model representing Material Stream"
  //1 -  Mixture, 2 - Liquid phase, 3 - Gas Phase
  
  
extends CrudeSimulator.GuessModels.InitialGuess ; 
  
  import CrudeSimulator.Files.*;
  parameter Integer Nc "Number of components" annotation(
    Dialog(tab = "Stream Specifications", group = "Component Parameters"));
  parameter CrudeSimulator.Files.CrudeDatabase.GeneralProperties C[Nc] "Component instances array" annotation(
    Dialog(tab = "Stream Specifications", group = "Component Parameters"));
  Real P(unit = "Pa", min = 0, start = Pg) "Pressure";
  Real T(unit = "K", start = Tg) "Temperature";
 Real xvap(unit = "-", start = xvapg, min = 0, max = 1) "Vapor Phase mole fraction";
  Real F_p[3](each unit = "mol/s", each min = 0, start={Fg,Fliqg,Fvapg}) "Total molar flow in phase";
  Real x_pc[3, Nc](each unit = "-", each min = 0, each max = 1, start={xguess,xg,yg}) "Component mole fraction in phase";
 Real Cp_p[3](each unit = "kJ/[kmol.K]",start={Hmixg,Hliqg,Hvapg}) "Phase molar specific heat";
  Real Cp_pc[3, Nc](each unit = "kJ/[kmol.K]") "Component molar specific heat in phase";
  Real H_p[3](each unit = "kJ/kmol",start={Hmixg,Hliqg,Hvapg}) "Phase molar enthalpy";
  Real H_pc[3, Nc](each unit = "kJ/kmol") "Component molar enthalpy in phase";
  CrudeSimulator.Files.Interfaces.matConn In(Nc = Nc) annotation(
    Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  CrudeSimulator.Files.Interfaces.matConn Out(Nc = Nc) annotation(
    Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));


equation
//Connector equations
  In.P = P;
  In.T = T;
  In.F = F_p[1];
  In.H = H_p[1];
  In.x_pc = x_pc;
  In.xvap = xvap;
  Out.P = P;
  Out.T = T;
  Out.F = F_p[1];
  Out.H = H_p[1];
  Out.x_pc = x_pc;
  Out.xvap = xvap;
//=====================================================================================
//Mole Balance
  F_p[1] = F_p[2] + F_p[3];
//Component Mole balance
 // x_pc[1, :] .* F_p[1] = x_pc[2, :] .* F_p[2] + x_pc[3, :] .* F_p[3];
//phase molar fraction
  xvap = F_p[3] / F_p[1];
//Energy Balance
  for i in 1:Nc loop
//Specific Heat and Enthalpy calculation
    Cp_pc[2, i] = Files.ThermodynamicFunctions.LiqCpId(C[i].LiqCp, T);
    Cp_pc[3, i] = Files.ThermodynamicFunctions.VapCpId(C[i].VapCp, T);
    H_pc[2, i] = Files.ThermodynamicFunctions.HLiqId(C[i].SH, C[i].VapCp, C[i].HOV, C[i].Tc, T);
    H_pc[3, i] = Files.ThermodynamicFunctions.HVapId(C[i].SH, C[i].VapCp, C[i].HOV, C[i].Tc, T);
  end for;
  for i in 2:3 loop
    Cp_p[i] = sum(x_pc[i, :] .* Cp_pc[i, :]) + Cpres_p[i];
    H_p[i] = sum(x_pc[i, :] .* H_pc[i, :]) + Hres_p[i];
  end for;
  Cp_p[1] = (1-xvap) * Cp_p[2] + xvap * Cp_p[3];
  Cp_pc[1, :] = x_pc[1, :] .* Cp_p[1];
  H_p[1] = (1-xvap) * H_p[2] + xvap * H_p[3];
  H_pc[1, :] = x_pc[1, :] .* H_p[1];
//VLE region
    for i in 1:Nc loop
      x_pc[3, i] = K_c[i] * x_pc[2, i];
      x_pc[2, i] = x_pc[1, i] ./ (1 + xvap * (K_c[i] - 1));
    end for;
    sum(x_pc[3, :]) = 1;
//sum y = 1

annotation(
    Documentation(info = "<html><head></head><body><div><!--StartFragment-->A <strong>Material Stream</strong> represents whatever enters and leaves the simulation passing through the unit operations.<!--EndFragment-->

</div><div><br></div><div>For variables which are decalared as 1-D array, the array size represent the phase where the array element indices 1 represents mixed phase, 2 represents liquid phase and 3 represents vapor phase.</div><div><br></div><div>For example, variable <b>F_p[3]</b> represents <i>Total molar flow in different phase</i>. So when simulated, the variables in the results will be as follow:</div><div>F_p[1] is Molar flow in mixed phase</div><div>F_p[2] is Molar flow in liquid phase</div><div>F_p[3] is Molar flow in vapor phase</div><div><br></div><div><br></div><div>For variables which are decalared as 2-D array, the first indice represent phase and second indice represents components.<div><br></div><div>For example, variable&nbsp;<b>F_pc[3,Nc]</b>&nbsp;represents <i>Component&nbsp;molar flow in different phase</i>. So when simulated, the variables in the results will be as follow:</div><div>F_pc[1,Nc] is Molar flow of Nc<sup>th</sup> in mixed phase</div><div>F_pc[2,Nc] is Molar flow of Nc<sup>th</sup> in liquid phase</div><div>F_pc[3,Nc] is Molar flow of Nc<sup>th</sup> in vapor phase</div></div><div><br></div><div><br></div><div><!--StartFragment--><span style=\"font-size: 12px;\">For demonstration on how to use this model to simulate a Material Stream,</span><span style=\"font-size: 12px;\">&nbsp;go to&nbsp;<a href=\"modelica://Simulator.Examples.CompositeMS\">Material Stream Example</a></span><!--EndFragment--></div></body></html>"));
     
end materialstream;
