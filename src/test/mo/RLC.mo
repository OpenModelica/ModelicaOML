model RLC
  Modelica.Electrical.Analog.Sources.SineVoltage sineVoltage(V = 1, freqHz = 50) annotation(Placement(transformation(extent = {{-10,-10},{10,10}}, rotation = 270, origin = {-40,22})));
  Modelica.Electrical.Analog.Basic.Resistor resistor(R = 1) annotation(Placement(transformation(extent = {{-8,50},{12,70}})));
  Modelica.Electrical.Analog.Basic.Capacitor capacitor(C = 0.00001) annotation(Placement(transformation(extent = {{-10,-10},{10,10}}, rotation = 270, origin = {52,-10})));
  Modelica.Electrical.Analog.Basic.Ground ground annotation(Placement(transformation(extent = {{-70,-26},{-50,-6}})));
  Modelica.Electrical.Analog.Basic.Inductor inductor(L = 0.01) annotation(
    Placement(visible = true, transformation(origin = {54, 28}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
equation
  connect(resistor.p, sineVoltage.p) annotation(
    Line(points = {{-8, 60}, {-40, 60}, {-40, 32}}, color = {0, 0, 255}, smooth = Smooth.None));
  connect(sineVoltage.n, capacitor.n) annotation(
    Line(points = {{-40, 12}, {-40, -20}, {52, -20}}, color = {0, 0, 255}, smooth = Smooth.None));
  connect(ground.p, sineVoltage.n) annotation(
    Line(points = {{-60, -6}, {-60, 12}, {-40, 12}}, color = {0, 0, 255}, smooth = Smooth.None));
  connect(inductor.n, capacitor.p) annotation(
    Line(points = {{54, 18}, {54, 9}, {52, 9}, {52, 0}}, color = {0, 0, 255}));
  connect(inductor.p, resistor.n) annotation(
    Line(points = {{54, 38}, {54, 60}, {12, 60}}, color = {0, 0, 255}));
  annotation(uses(Modelica(version = "3.2")), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}), graphics), experiment(StopTime = 0.1), experimentSetupOutput);
end RLC;
