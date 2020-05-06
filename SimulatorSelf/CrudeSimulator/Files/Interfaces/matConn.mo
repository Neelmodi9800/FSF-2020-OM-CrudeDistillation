within CrudeSimulator.Files.Interfaces;

connector matConn
 Real P, T, F, H, x_pc[3, Nc], xvap, S;
    parameter Integer Nc;
  annotation(
    Icon(coordinateSystem(initialScale = 0.1), graphics = {Rectangle(lineColor = {0, 70, 70}, fillColor = {0, 70, 70}, fillPattern = FillPattern.Solid, extent = {{-50, 50}, {50, -50}})}));
end matConn;
