package trial
  model matstream
    extends CrudeSimulator.Streams.materialstream;
    extends CrudeSimulator.Files.ThermodynamicPackages.RaoultsLaw;
  end matstream;

  model check
    parameter Integer Nc = 4;
    import data = CrudeSimulator.Files.CrudeDatabase;
    parameter data.Pseudocomponent1 p1;
    parameter data.Pseudocomponent2 p2;
    parameter data.Pseudocomponent3 p3;
    parameter data.Pseudocomponent4 p4;
    parameter CrudeSimulator.Files.CrudeDatabase.GeneralProperties C[Nc] = {p1, p2, p3, p4};
    parameter Real Tcheck = 400;
    matstream S01(Nc = Nc, C = C) annotation(
      Placement(visible = true, transformation(origin = {-60, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    S01.F_p[1] = 50;
    S01.T = Tcheck;
    S01.x_pc[1, :] = {0.25, 0.2, 0.3, 0.25};
    S01.P = 101325;
  end check;
end trial;
