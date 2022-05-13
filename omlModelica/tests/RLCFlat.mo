class RLC
  parameter Real sineVoltage.V(quantity = "ElectricPotential", unit = "V", start = 1.0) = 1.0 "Amplitude of sine wave";
  parameter Real sineVoltage.phase(quantity = "Angle", unit = "rad", displayUnit = "deg") = 0.0 "Phase of sine wave";
  parameter Real sineVoltage.freqHz(quantity = "Frequency", unit = "Hz", start = 1.0) = 50.0 "Frequency of sine wave";
  Real sineVoltage.v(quantity = "ElectricPotential", unit = "V") "Voltage drop of the two pins (= p.v - n.v)";
  Real sineVoltage.i(quantity = "ElectricCurrent", unit = "A") "Current flowing from pin p to pin n";
  Real sineVoltage.p.v(quantity = "ElectricPotential", unit = "V") "Potential at the pin";
  Real sineVoltage.p.i(quantity = "ElectricCurrent", unit = "A") "Current flowing into the pin";
  Real sineVoltage.n.v(quantity = "ElectricPotential", unit = "V") "Potential at the pin";
  Real sineVoltage.n.i(quantity = "ElectricCurrent", unit = "A") "Current flowing into the pin";
  parameter Real sineVoltage.offset(quantity = "ElectricPotential", unit = "V") = 0.0 "Voltage offset";
  parameter Real sineVoltage.startTime(quantity = "Time", unit = "s") = 0.0 "Time offset";
  final parameter Real sineVoltage.signalSource.amplitude = sineVoltage.V "Amplitude of sine wave";
  final parameter Real sineVoltage.signalSource.freqHz(quantity = "Frequency", unit = "Hz", start = 1.0) = sineVoltage.freqHz "Frequency of sine wave";
  final parameter Real sineVoltage.signalSource.phase(quantity = "Angle", unit = "rad", displayUnit = "deg") = sineVoltage.phase "Phase of sine wave";
  Real sineVoltage.signalSource.y "Connector of Real output signal";
  final parameter Real sineVoltage.signalSource.offset = sineVoltage.offset "Offset of output signal y";
  final parameter Real sineVoltage.signalSource.startTime(quantity = "Time", unit = "s") = sineVoltage.startTime "Output y = offset for time < startTime";
  
  parameter Real resistor.R(quantity = "Resistance", unit = "Ohm", start = 1.0) = 1.0 "Resistance at temperature T_ref";
  parameter Real resistor.T_ref(quantity = "ThermodynamicTemperature", unit = "K", displayUnit = "degC", min = 0.0, start = 288.15, nominal = 300.0) = 300.15 "Reference temperature";
  parameter Real resistor.alpha(quantity = "LinearTemperatureCoefficient", unit = "1/K") = 0.0 "Temperature coefficient of resistance (R_actual = R*(1 + alpha*(T_heatPort - T_ref))";
  Real resistor.v(quantity = "ElectricPotential", unit = "V") "Voltage drop of the two pins (= p.v - n.v)";
  Real resistor.i(quantity = "ElectricCurrent", unit = "A") "Current flowing from pin p to pin n";
  Real resistor.p.v(quantity = "ElectricPotential", unit = "V") "Potential at the pin";
  Real resistor.p.i(quantity = "ElectricCurrent", unit = "A") "Current flowing into the pin";
  Real resistor.n.v(quantity = "ElectricPotential", unit = "V") "Potential at the pin";
  Real resistor.n.i(quantity = "ElectricCurrent", unit = "A") "Current flowing into the pin";
  final parameter Boolean resistor.useHeatPort = false "=true, if heatPort is enabled";
  parameter Real resistor.T(quantity = "ThermodynamicTemperature", unit = "K", displayUnit = "degC", min = 0.0, start = 288.15, nominal = 300.0) = resistor.T_ref "Fixed device temperature if useHeatPort = false";
  Real resistor.LossPower(quantity = "Power", unit = "W") "Loss power leaving component via heatPort";
  Real resistor.T_heatPort(quantity = "ThermodynamicTemperature", unit = "K", displayUnit = "degC", min = 0.0, start = 288.15, nominal = 300.0) "Temperature of heatPort";
  Real resistor.R_actual(quantity = "Resistance", unit = "Ohm") "Actual resistance = R*(1 + alpha*(T_heatPort - T_ref))";
  Real capacitor.v(quantity = "ElectricPotential", unit = "V", start = 0.0) "Voltage drop of the two pins (= p.v - n.v)";
  Real capacitor.i(quantity = "ElectricCurrent", unit = "A") "Current flowing from pin p to pin n";
  Real capacitor.p.v(quantity = "ElectricPotential", unit = "V") "Potential at the pin";
  Real capacitor.p.i(quantity = "ElectricCurrent", unit = "A") "Current flowing into the pin";
  Real capacitor.n.v(quantity = "ElectricPotential", unit = "V") "Potential at the pin";
  Real capacitor.n.i(quantity = "ElectricCurrent", unit = "A") "Current flowing into the pin";
  parameter Real capacitor.C(quantity = "Capacitance", unit = "F", min = 0.0, start = 1.0) = 1e-05 "Capacitance";
  Real inductor.v(quantity = "ElectricPotential", unit = "V") "Voltage drop of the two pins (= p.v - n.v)";
  Real inductor.i(quantity = "ElectricCurrent", unit = "A", start = 0.0) "Current flowing from pin p to pin n";
  Real inductor.p.v(quantity = "ElectricPotential", unit = "V") "Potential at the pin";
  Real inductor.p.i(quantity = "ElectricCurrent", unit = "A") "Current flowing into the pin";
  Real inductor.n.v(quantity = "ElectricPotential", unit = "V") "Potential at the pin";
  Real inductor.n.i(quantity = "ElectricCurrent", unit = "A") "Current flowing into the pin";
  parameter Real inductor.L(quantity = "Inductance", unit = "H", start = 1.0) = 0.01 "Inductance";
  Real ground.p.v(quantity = "ElectricPotential", unit = "V") "Potential at the pin";
  Real ground.p.i(quantity = "ElectricCurrent", unit = "A") "Current flowing into the pin";
equation
  inductor.p.v = resistor.n.v;
  resistor.p.v = sineVoltage.p.v;
  inductor.n.v = capacitor.p.v;
  ground.p.v = sineVoltage.n.v;
  ground.p.v = capacitor.n.v;
  resistor.p.i + sineVoltage.p.i = 0.0;
  ground.p.i + capacitor.n.i + sineVoltage.n.i = 0.0;
  inductor.p.i + resistor.n.i = 0.0;
  inductor.n.i + capacitor.p.i = 0.0;
  sineVoltage.signalSource.y = sineVoltage.signalSource.offset + (if time < sineVoltage.signalSource.startTime then 0.0 else sineVoltage.signalSource.amplitude * sin(6.283185307179586 * sineVoltage.signalSource.freqHz * (time - sineVoltage.signalSource.startTime) + sineVoltage.signalSource.phase));
  sineVoltage.v = sineVoltage.signalSource.y;
  sineVoltage.v = sineVoltage.p.v - sineVoltage.n.v;
  0.0 = sineVoltage.p.i + sineVoltage.n.i;
  sineVoltage.i = sineVoltage.p.i;
  assert(1.0 + resistor.alpha * (resistor.T_heatPort - resistor.T_ref) >= 1e-15, "Temperature outside scope of model!");
  resistor.R_actual = resistor.R * (1.0 + resistor.alpha * (resistor.T_heatPort - resistor.T_ref));
  resistor.v = resistor.R_actual * resistor.i;
  resistor.LossPower = resistor.v * resistor.i;
  resistor.T_heatPort = resistor.T;
  resistor.v = resistor.p.v - resistor.n.v;
  0.0 = resistor.p.i + resistor.n.i;
  resistor.i = resistor.p.i;
  capacitor.i = capacitor.C * der(capacitor.v);
  capacitor.v = capacitor.p.v - capacitor.n.v;
  0.0 = capacitor.p.i + capacitor.n.i;
  capacitor.i = capacitor.p.i;
  inductor.L * der(inductor.i) = inductor.v;
  inductor.v = inductor.p.v - inductor.n.v;
  0.0 = inductor.p.i + inductor.n.i;
  inductor.i = inductor.p.i;
  ground.p.v = 0.0;
end RLC;
