within Distillation;

model matstream_dummy

  parameter Integer Nc "Number of components" annotation(
    Dialog(tab = "Stream Specifications", group = "Component Parameters"));
  Real F_p(unit = "mol/s", min = 0) "Total molar flow in mixure";
  Real x_pc[Nc](each unit = " - ", each min = 0, each max = 1) "Mole fraction of component in mixure";
  matconn_sim In(Nc = Nc) annotation(
    Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  matconn_sim Out(Nc = Nc) annotation(
    Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
 
equation
  
  In.F = F_p;
  In.x_pc = x_pc;
  Out.F = F_p;
  Out.x_pc = x_pc;

 
end matstream_dummy;
