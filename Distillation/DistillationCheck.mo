within Distillation;

package DistillationCheck
  model matstream
    extends matstream_dummy;
  end matstream;

      model flowsheet
              parameter Integer Nc = 4;
              Distillation.DistillationColumn dis1(Nc = Nc, NI = 6, In_s = 3) annotation(
                      Placement(visible = true, transformation(origin = {2, 16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
              matstream S01(Nc = Nc) annotation(
                      Placement(visible = true, transformation(origin = {-74, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
              matstream S03(Nc = Nc) annotation(
                      Placement(visible = true, transformation(origin = {74, -22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
              matstream S02(Nc = Nc) annotation(
                      Placement(visible = true, transformation(origin = {68, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
              connect(dis1.In, S01.Out) annotation(
                      Line(points = {{-22, 16}, {-54, 16}, {-54, 20}, {-64, 20}, {-64, 20}}, color = {0, 70, 70}));
              connect(dis1.Bot, S03.In) annotation(
                      Line(points = {{27, -14}, {51.5, -14}, {51.5, -22}, {64, -22}}, color = {0, 70, 70}));
              connect(dis1.Dist, S02.In) annotation(
                      Line(points = {{27, 46}, {58, 46}}, color = {0, 70, 70}));
              S01.F_p = 10;
              S01.x_pc = {0.25, 0.72, 0.03, 0.02};
              dis1.k_c[:] = {0.0852, 0.8969, 9.309, 34.42};
              dis1.RR = 2;
              S03.F_p = 3;
      end flowsheet;
  
  
end DistillationCheck;
