package ThermoPower "Open library for thermal power plant simulation"
  extends Modelica.Icons.Package;
  import SI = Modelica.SIunits;
  import NonSI = Modelica.SIunits.Conversions.NonSIunits;

  model System "System wide properties and defaults"
    parameter Boolean allowFlowReversal = true "= false to restrict to design flow direction (flangeA -> flangeB)" annotation(Evaluate = true);
    parameter Choices.Init.Options initOpt = ThermoPower.Choices.Init.Options.fixedState;
    parameter SI.Pressure p_amb = 101325 "Ambient pressure";
    parameter SI.Temperature T_amb = 293.15 "Ambient Temperature (dry bulb)";
    parameter SI.Temperature T_wb = 288.15 "Ambient temperature (wet bulb)";
    parameter SI.Frequency fnom = 50 "Nominal grid frequency";
    annotation(defaultComponentPrefixes = "inner", missingInnerMessage = "The System object is missing, please drag it on the top layer of your model");
  end System;

  package Gas "Models of components with ideal gases as working fluid"
    connector Flange "Flange connector for gas flows"
      replaceable package Medium = Modelica.Media.Interfaces.PartialMedium;
      flow Medium.MassFlowRate m_flow "Mass flow rate from the connection point into the component";
      Medium.AbsolutePressure p "Thermodynamic pressure in the connection point";
      stream Medium.SpecificEnthalpy h_outflow "Specific thermodynamic enthalpy close to the connection point if m_flow < 0";
      stream Medium.MassFraction[Medium.nXi] Xi_outflow "Independent mixture mass fractions m_i/m close to the connection point if m_flow < 0";
      stream Medium.ExtraProperty[Medium.nC] C_outflow "Properties c_i/m close to the connection point if m_flow < 0";
    end Flange;

    connector FlangeA "A-type flange connector for gas flows"
      extends Flange;
    end FlangeA;

    connector FlangeB "B-type flange connector for gas flows"
      extends Flange;
    end FlangeB;

    extends Modelica.Icons.Package;

    model SourcePressure "Pressure source for gas flows"
      extends Icons.Gas.SourceP;
      replaceable package Medium = Modelica.Media.Interfaces.PartialMedium annotation(choicesAllMatching = true);
      Medium.BaseProperties gas(p(start = p0), T(start = T), Xi(start = Xnom[1:Medium.nXi]));
      parameter SI.Pressure p0 = 101325 "Nominal pressure";
      parameter Units.HydraulicResistance R = 0 "Hydraulic resistance";
      parameter Medium.Temperature T = 300 "Nominal temperature";
      parameter Medium.MassFraction[Medium.nX] Xnom = Medium.reference_X "Nominal gas composition";
      parameter Boolean allowFlowReversal = system.allowFlowReversal "= true to allow flow reversal, false restricts to design direction" annotation(Evaluate = true);
      parameter Boolean use_in_p0 = false "Use connector input for the pressure";
      parameter Boolean use_in_T = false "Use connector input for the temperature";
      parameter Boolean use_in_X = false "Use connector input for the composition";
      outer ThermoPower.System system "System wide properties";
      FlangeB flange(redeclare package Medium = Medium, m_flow(max = if allowFlowReversal then +Modelica.Constants.inf else 0));
      Modelica.Blocks.Interfaces.RealInput in_p0 if use_in_p0;
      Modelica.Blocks.Interfaces.RealInput in_T if use_in_T;
      Modelica.Blocks.Interfaces.RealInput[Medium.nX] in_X if use_in_X;
    protected
      Modelica.Blocks.Interfaces.RealInput in_p0_internal;
      Modelica.Blocks.Interfaces.RealInput in_T_internal;
      Modelica.Blocks.Interfaces.RealInput[Medium.nX] in_X_internal;
    equation
      if R > 0 then
        flange.p = gas.p + flange.m_flow * R;
      else
        flange.p = gas.p;
      end if;
      gas.p = in_p0_internal;
      if not use_in_p0 then
        in_p0_internal = p0 "Pressure set by parameter";
      end if;
      gas.T = in_T_internal;
      if not use_in_T then
        in_T_internal = T "Temperature set by parameter";
      end if;
      gas.Xi = in_X_internal[1:Medium.nXi];
      if not use_in_X then
        in_X_internal = Xnom "Composition set by parameter";
      end if;
      flange.h_outflow = gas.h;
      flange.Xi_outflow = gas.Xi;
      connect(in_p0, in_p0_internal);
      connect(in_T, in_T_internal);
      connect(in_X, in_X_internal);
    end SourcePressure;

    model SinkPressure "Pressure sink for gas flows"
      extends Icons.Gas.SourceP;
      replaceable package Medium = Modelica.Media.Interfaces.PartialMedium annotation(choicesAllMatching = true);
      Medium.BaseProperties gas(p(start = p0), T(start = T), Xi(start = Xnom[1:Medium.nXi]));
      parameter Medium.AbsolutePressure p0 = 101325 "Nominal pressure";
      parameter Medium.Temperature T = 300 "Nominal temperature";
      parameter Medium.MassFraction[Medium.nX] Xnom = Medium.reference_X "Nominal gas composition";
      parameter Units.HydraulicResistance R = 0 "Hydraulic Resistance";
      parameter Boolean allowFlowReversal = system.allowFlowReversal "= true to allow flow reversal, false restricts to design direction" annotation(Evaluate = true);
      parameter Boolean use_in_p0 = false "Use connector input for the pressure";
      parameter Boolean use_in_T = false "Use connector input for the temperature";
      parameter Boolean use_in_X = false "Use connector input for the composition";
      outer ThermoPower.System system "System wide properties";
      FlangeA flange(redeclare package Medium = Medium, m_flow(min = if allowFlowReversal then -Modelica.Constants.inf else 0));
      Modelica.Blocks.Interfaces.RealInput in_p0 if use_in_p0;
      Modelica.Blocks.Interfaces.RealInput in_T if use_in_T;
      Modelica.Blocks.Interfaces.RealInput[Medium.nX] in_X if use_in_X;
    protected
      Modelica.Blocks.Interfaces.RealInput in_p0_internal;
      Modelica.Blocks.Interfaces.RealInput in_T_internal;
      Modelica.Blocks.Interfaces.RealInput[Medium.nX] in_X_internal;
    equation
      if R > 0 then
        flange.p = gas.p + flange.m_flow * R;
      else
        flange.p = gas.p;
      end if;
      gas.p = in_p0_internal;
      if not use_in_p0 then
        in_p0_internal = p0 "Pressure set by parameter";
      end if;
      gas.T = in_T_internal;
      if not use_in_T then
        in_T_internal = T "Temperature set by parameter";
      end if;
      gas.Xi = in_X_internal[1:Medium.nXi];
      if not use_in_X then
        in_X_internal = Xnom "Composition set by parameter";
      end if;
      flange.h_outflow = gas.h;
      flange.Xi_outflow = gas.Xi;
      connect(in_p0, in_p0_internal);
      connect(in_T, in_T_internal);
      connect(in_X, in_X_internal);
    end SinkPressure;

    model SourceMassFlow "Flow rate source for gas flows"
      extends Icons.Gas.SourceW;
      replaceable package Medium = Modelica.Media.Interfaces.PartialMedium annotation(choicesAllMatching = true);
      Medium.BaseProperties gas(p(start = p0), T(start = T), Xi(start = Xnom[1:Medium.nXi]));
      parameter Medium.AbsolutePressure p0 = 101325 "Nominal pressure";
      parameter Medium.Temperature T = 300 "Nominal temperature";
      parameter Medium.MassFraction[Medium.nX] Xnom = Medium.reference_X "Nominal gas composition";
      parameter Medium.MassFlowRate w0 = 0 "Nominal mass flowrate";
      parameter Units.HydraulicConductance G = 0 "HydraulicConductance";
      parameter Boolean allowFlowReversal = system.allowFlowReversal "= true to allow flow reversal, false restricts to design direction" annotation(Evaluate = true);
      parameter Boolean use_in_w0 = false "Use connector input for the nominal flow rate";
      parameter Boolean use_in_T = false "Use connector input for the temperature";
      parameter Boolean use_in_X = false "Use connector input for the composition";
      outer ThermoPower.System system "System wide properties";
      Medium.MassFlowRate w "Nominal mass flow rate";
      FlangeB flange(redeclare package Medium = Medium, m_flow(max = if allowFlowReversal then +Modelica.Constants.inf else 0));
      Modelica.Blocks.Interfaces.RealInput in_w0 if use_in_w0;
      Modelica.Blocks.Interfaces.RealInput in_T if use_in_T;
      Modelica.Blocks.Interfaces.RealInput[Medium.nX] in_X if use_in_X;
    protected
      Modelica.Blocks.Interfaces.RealInput in_w0_internal;
      Modelica.Blocks.Interfaces.RealInput in_T_internal;
      Modelica.Blocks.Interfaces.RealInput[Medium.nX] in_X_internal;
    equation
      if G > 0 then
        flange.m_flow = (-w) + (flange.p - p0) * G;
      else
        flange.m_flow = -w;
      end if;
      w = in_w0_internal;
      if not use_in_w0 then
        in_w0_internal = w0 "Flow rate set by parameter";
      end if;
      gas.T = in_T_internal;
      if not use_in_T then
        in_T_internal = T "Temperature set by parameter";
      end if;
      gas.Xi = in_X_internal[1:Medium.nXi];
      if not use_in_X then
        in_X_internal = Xnom "Composition set by parameter";
      end if;
      flange.p = gas.p;
      flange.h_outflow = gas.h;
      flange.Xi_outflow = gas.Xi;
      connect(in_w0, in_w0_internal);
      connect(in_T, in_T_internal);
      connect(in_X, in_X_internal);
    end SourceMassFlow;

    model Flow1DFV "1-dimensional fluid flow model for gas (finite volumes)"
      extends BaseClasses.Flow1DBase;
      Thermal.DHTVolumes wall(final N = Nw);
      replaceable model HeatTransfer = Thermal.HeatTransferFV.IdealHeatTransfer constrainedby ThermoPower.Thermal.BaseClasses.DistributedHeatTransferFV;
      HeatTransfer heatTransfer(redeclare package Medium = Medium, final Nf = N, final Nw = Nw, final Nt = Nt, final L = L, final A = A, final Dhyd = Dhyd, final omega = omega, final wnom = wnom / Nt, final w = w * ones(N), final fluidState = gas.state) "Instantiated heat transfer model";
      parameter SI.PerUnit wnm = 1e-2 "Maximum fraction of the nominal flow rate allowed as reverse flow";
      parameter Boolean fixedMassFlowSimplified = false "Fix flow rate = wnom for simplified homotopy model";
      Medium.BaseProperties[N] gas "Gas nodal properties";
      SI.Pressure Dpfric "Pressure drop due to friction";
      SI.Length omega_hyd "Wet perimeter (single tube)";
      Real Kf "Friction factor";
      Real Kfl "Linear friction factor";
      Real dwdt "Time derivative of mass flow rate";
      SI.PerUnit Cf "Fanning friction factor";
      Medium.MassFlowRate w(start = wnom / Nt) "Mass flowrate (single tube)";
      SI.Temperature[N - 1] Ttilde(start = Tstart[2:N], each stateSelect = StateSelect.prefer) "Temperature state variables";
      Medium.Temperature[N] T(start = Tstart) "Node temperatures";
      Medium.SpecificEnthalpy[N] h "Node specific enthalpies";
      Medium.Temperature Tin(start = Tstartin);
      Medium.MassFraction[if UniformComposition or Medium.fixedX then 1 else N - 1, nX] Xtilde(start = ones(size(Xtilde, 1), size(Xtilde, 2)) * diagonal(Xstart[1:nX]), each stateSelect = StateSelect.prefer) "Composition state variables";
      Medium.MassFlowRate[N - 1] wbar(each start = wnom / Nt);
      SI.Power[N - 1] Q_single = heatTransfer.Qvol / Nt "Heat flows entering the volumes from the lateral boundary (single tube)";
      SI.Velocity[N] u "Fluid velocity";
      Medium.AbsolutePressure p(start = pstart, stateSelect = StateSelect.prefer);
      SI.Time Tr "Residence time";
      SI.Mass M "Gas Mass (single tube)";
      SI.Mass Mtot "Gas Mass (total)";
      SI.Power Q "Total heat flow through the wall (all Nt tubes)";
    protected
      parameter SI.Length l = L / (N - 1) "Length of a single volume";
      Medium.Density[N - 1] rhobar "Fluid average density";
      SI.SpecificVolume[N - 1] vbar "Fluid average specific volume";
      Medium.DerDensityByPressure[N - 1] drbdp "Derivative of average density by pressure";
      Medium.DerDensityByTemperature[N - 1] drbdT1 "Derivative of average density by left temperature";
      Medium.DerDensityByTemperature[N - 1] drbdT2 "Derivative of average density by right temperature";
      Real[N - 1, nX] drbdX1(each unit = "kg/m3") "Derivative of average density by left composition";
      Real[N - 1, nX] drbdX2(each unit = "kg/m3") "Derivative of average density by right composition";
      Medium.SpecificHeatCapacity[N - 1] cvbar "Average cv";
      SI.MassFlowRate[N - 1] dMdt "Derivative of mass in a finite volume";
      Medium.SpecificHeatCapacity[N] cv;
      Medium.DerDensityByTemperature[N] dddT "Derivative of density by temperature";
      Medium.DerDensityByPressure[N] dddp "Derivative of density by pressure";
      Real[N, nX] dddX(each unit = "kg/m3") "Derivative of density by composition";
    initial equation
      if initOpt == Choices.Init.Options.noInit or QuasiStatic then
      elseif initOpt == Choices.Init.Options.fixedState then
        if not noInitialPressure then
          p = pstart;
        end if;
        Ttilde = Tstart[2:N];
        if not Medium.fixedX then
          Xtilde = ones(size(Xtilde, 1), size(Xtilde, 2)) * diagonal(Xstart[1:nX]);
        end if;
      elseif initOpt == Choices.Init.Options.steadyState then
        if not Medium.singleState and not noInitialPressure then
          der(p) = 0;
        end if;
        der(Ttilde) = zeros(N - 1);
        if not Medium.fixedX then
          der(Xtilde) = zeros(size(Xtilde, 1), size(Xtilde, 2));
        end if;
      elseif initOpt == Choices.Init.Options.steadyStateNoP then
        der(Ttilde) = zeros(N - 1);
        if not Medium.fixedX then
          der(Xtilde) = zeros(size(Xtilde, 1), size(Xtilde, 2));
        end if;
      else
        assert(false, "Unsupported initialisation option");
      end if;
    equation
      assert(FFtype == ThermoPower.Choices.Flow1D.FFtypes.NoFriction or dpnom > 0, "dpnom=0 not supported, it is also used in the homotopy trasformation during the inizialization");
      omega_hyd = 4 * A / Dhyd;
      if FFtype == ThermoPower.Choices.Flow1D.FFtypes.Kfnom then
        Kf = Kfnom * Kfc;
        Cf = 2 * Kf * A ^ 3 / (omega_hyd * L);
      elseif FFtype == ThermoPower.Choices.Flow1D.FFtypes.OpPoint then
        Kf = dpnom * rhonom / (wnom / Nt) ^ 2 * Kfc;
        Cf = 2 * Kf * A ^ 3 / (omega_hyd * L);
      elseif FFtype == ThermoPower.Choices.Flow1D.FFtypes.Cfnom then
        Kf = Cfnom * omega_hyd * L / (2 * A ^ 3) * Kfc;
        Cf = Cfnom * Kfc;
      elseif FFtype == ThermoPower.Choices.Flow1D.FFtypes.Colebrook then
        Cf = f_colebrook(w, Dhyd / A, e, Medium.dynamicViscosity(gas[integer(N / 2)].state)) * Kfc;
        Kf = Cf * omega_hyd * L / (2 * A ^ 3);
      elseif FFtype == ThermoPower.Choices.Flow1D.FFtypes.NoFriction then
        Cf = 0;
        Kf = 0;
      else
        assert(false, "Unsupported FFtype");
        Cf = 0;
        Kf = 0;
      end if;
      assert(Kf >= 0, "Negative friction coefficient");
      Kfl = wnom / Nt * wnf * Kf "Linear friction factor";
      dwdt = if DynamicMomentum and not QuasiStatic then der(w) else 0;
      sum(dMdt) = (infl.m_flow + outfl.m_flow) / Nt "Mass balance";
      L / A * dwdt + outfl.p - infl.p + Dpfric = 0 "Momentum balance";
      Dpfric = if FFtype == ThermoPower.Choices.Flow1D.FFtypes.NoFriction then 0 else homotopy(smooth(1, Kf * squareReg(w, wnom / Nt * wnf)) * sum(vbar) / (N - 1), dpnom / (wnom / Nt) * w) "Pressure drop due to friction";
      for j in 1:N - 1 loop
        if not QuasiStatic then
          A * l * rhobar[j] * cvbar[j] * der(Ttilde[j]) + wbar[j] * (gas[j + 1].h - gas[j].h) = Q_single[j] "Energy balance";
          dMdt[j] = A * l * (drbdp[j] * der(p) + drbdT1[j] * der(gas[j].T) + drbdT2[j] * der(gas[j + 1].T) + vector(drbdX1[j, :]) * vector(der(gas[j].X)) + vector(drbdX2[j, :]) * vector(der(gas[j + 1].X))) "Mass balance";
          if avoidInletEnthalpyDerivative and j == 1 then
            rhobar[j] = gas[j + 1].d;
            drbdp[j] = dddp[j + 1];
            drbdT1[j] = 0;
            drbdT2[j] = dddT[j + 1];
            drbdX1[j, :] = zeros(size(Xtilde, 2));
            drbdX2[j, :] = dddX[j + 1, :];
          else
            rhobar[j] = (gas[j].d + gas[j + 1].d) / 2;
            drbdp[j] = (dddp[j] + dddp[j + 1]) / 2;
            drbdT1[j] = dddT[j] / 2;
            drbdT2[j] = dddT[j + 1] / 2;
            drbdX1[j, :] = dddX[j, :] / 2;
            drbdX2[j, :] = dddX[j + 1, :] / 2;
          end if;
          vbar[j] = 1 / rhobar[j];
          if fixedMassFlowSimplified then
            wbar[j] = homotopy(infl.m_flow / Nt - sum(dMdt[1:j - 1]) - dMdt[j] / 2, wnom / Nt);
          else
            wbar[j] = infl.m_flow / Nt - sum(dMdt[1:j - 1]) - dMdt[j] / 2;
          end if;
          cvbar[j] = (cv[j] + cv[j + 1]) / 2;
        else
          wbar[j] * (gas[j + 1].h - gas[j].h) = Q_single[j] "Energy balance";
          dMdt[j] = 0 "Mass balance";
          rhobar[j] = 0;
          drbdp[j] = 0;
          drbdT1[j] = 0;
          drbdT2[j] = 0;
          drbdX1[j, :] = zeros(nX);
          drbdX2[j, :] = zeros(nX);
          vbar[j] = 0;
          if fixedMassFlowSimplified then
            wbar[j] = homotopy(infl.m_flow / Nt, wnom / Nt);
          else
            wbar[j] = infl.m_flow / Nt;
          end if;
          cvbar[j] = 0;
        end if;
      end for;
      Q = heatTransfer.Q "Total heat flow through the lateral boundary";
      if Medium.fixedX then
        Xtilde = fill(Medium.reference_X, 1);
      elseif QuasiStatic then
        Xtilde = fill(gas[1].X, size(Xtilde, 1)) "Gas composition equal to actual inlet";
      elseif UniformComposition then
        der(Xtilde[1, :]) = homotopy(1 / L * sum(u) / N * (gas[1].X - gas[N].X), 1 / L * unom * (gas[1].X - gas[N].X)) "Partial mass balance for the whole pipe";
      else
        for j in 1:N - 1 loop
          der(Xtilde[j, :]) = homotopy((u[j + 1] + u[j]) / (2 * l) * (gas[j].X - gas[j + 1].X), 1 / L * unom * (gas[j].X - gas[j + 1].X)) "Partial mass balance for single volume";
        end for;
      end if;
      for j in 1:N loop
        u[j] = w / (gas[j].d * A) "Gas velocity";
        gas[j].p = p;
        gas[j].T = T[j];
        gas[j].h = h[j];
      end for;
      for j in 1:N loop
        if not QuasiStatic then
          cv[j] = Medium.heatCapacity_cv(gas[j].state);
          dddT[j] = Medium.density_derT_p(gas[j].state);
          dddp[j] = Medium.density_derp_T(gas[j].state);
          if nX > 0 then
            dddX[j, :] = Medium.density_derX(gas[j].state);
          end if;
        else
          cv[j] = 0;
          dddT[j] = 0;
          dddp[j] = 0;
          dddX[j, :] = zeros(nX);
        end if;
      end for;
      if HydraulicCapacitance == ThermoPower.Choices.Flow1D.HCtypes.Upstream then
        p = infl.p;
        w = -outfl.m_flow / Nt;
      else
        p = outfl.p;
        w = infl.m_flow / Nt;
      end if;
      infl.h_outflow = gas[1].h;
      outfl.h_outflow = gas[N].h;
      infl.Xi_outflow = gas[1].Xi;
      outfl.Xi_outflow = gas[N].Xi;
      gas[1].h = inStream(infl.h_outflow);
      gas[2:N].T = Ttilde;
      gas[1].Xi = inStream(infl.Xi_outflow);
      for j in 2:N loop
        gas[j].Xi = Xtilde[if UniformComposition then 1 else j - 1, 1:nXi];
      end for;
      connect(wall, heatTransfer.wall);
      Tin = gas[1].T;
      M = sum(rhobar) * A * l "Fluid mass (single tube)";
      Mtot = M * Nt "Fluid mass (total)";
      Tr = noEvent(M / max(infl.m_flow / Nt, Modelica.Constants.eps)) "Residence time";
      assert(infl.m_flow > (-wnom * wnm), "Reverse flow not allowed, maybe you connected the component with wrong orientation");
    end Flow1DFV;

    function f_colebrook "Fanning friction factor for water/steam flows"
      input SI.MassFlowRate w;
      input Real D_A;
      input Real e;
      input SI.DynamicViscosity mu;
      output SI.PerUnit f;
    protected
      Real Re;
    algorithm
      Re := w * D_A / mu;
      Re := if Re > 2100 then Re else 2100;
      f := 0.332 / log(e / 3.7 + 5.47 / Re ^ 0.9) ^ 2;
    end f_colebrook;

    package BaseClasses
      extends Modelica.Icons.BasesPackage;

      partial model Flow1DBase "Basic interface for 1-dimensional water/steam fluid flow models"
        extends Icons.Gas.Tube;
        import ThermoPower.Choices.Flow1D.FFtypes;
        import ThermoPower.Choices.Flow1D.HCtypes;
        replaceable package Medium = Modelica.Media.Interfaces.PartialMedium annotation(choicesAllMatching = true);
        parameter Integer N(min = 2) = 2 "Number of nodes for thermal variables";
        parameter Integer Nw = N - 1 "Number of volumes on the wall interface";
        parameter Integer Nt = 1 "Number of tubes in parallel";
        parameter SI.Distance L "Tube length";
        parameter SI.Position H = 0 "Elevation of outlet over inlet";
        parameter SI.Area A "Cross-sectional area (single tube)";
        parameter SI.Length omega "Perimeter of heat transfer surface (single tube)";
        parameter SI.Length Dhyd "Hydraulic Diameter (single tube)";
        parameter Medium.MassFlowRate wnom "Nominal mass flowrate (total)";
        parameter ThermoPower.Choices.Flow1D.FFtypes FFtype = ThermoPower.Choices.Flow1D.FFtypes.NoFriction "Friction Factor Type" annotation(Evaluate = true);
        parameter SI.PressureDifference dpnom = 0 "Nominal pressure drop";
        parameter Real Kfnom = 0 "Nominal hydraulic resistance coefficient";
        parameter Medium.Density rhonom = 0 "Nominal inlet density";
        parameter SI.PerUnit Cfnom = 0 "Nominal Fanning friction factor";
        parameter SI.PerUnit e = 0 "Relative roughness (ratio roughness/diameter)";
        parameter Real Kfc = 1 "Friction factor correction coefficient";
        parameter Boolean DynamicMomentum = false "Inertial phenomena accounted for" annotation(Evaluate = true);
        parameter Boolean UniformComposition = true "Uniform gas composition is assumed" annotation(Evaluate = true);
        parameter Boolean QuasiStatic = false "Quasi-static model (mass, energy and momentum static balances" annotation(Evaluate = true);
        parameter HCtypes HydraulicCapacitance = HCtypes.Downstream "1: Upstream, 2: Downstream";
        parameter Boolean avoidInletEnthalpyDerivative = true "Avoid inlet enthalpy derivative";
        parameter Boolean allowFlowReversal = system.allowFlowReversal "= true to allow flow reversal, false restricts to design direction" annotation(Evaluate = true);
        outer ThermoPower.System system "System wide properties";
        parameter Medium.AbsolutePressure pstart = 1e5 "Pressure start value";
        parameter Medium.Temperature Tstartbar = 300 "Avarage temperature start value";
        parameter Medium.Temperature Tstartin = Tstartbar "Inlet temperature start value";
        parameter Medium.Temperature Tstartout = Tstartbar "Outlet temperature start value";
        parameter Medium.Temperature[N] Tstart = linspace(Tstartin, Tstartout, N) "Start value of temperature vector (initialized by default)";
        final parameter SI.Velocity unom = 10 "Nominal velocity for simplified equation";
        parameter Real wnf = 0.01 "Fraction of nominal flow rate at which linear friction equals turbulent friction";
        parameter Medium.MassFraction[nX] Xstart = Medium.reference_X "Start gas composition";
        parameter Choices.Init.Options initOpt = system.initOpt "Initialisation option";
        parameter Boolean noInitialPressure = false "Remove initial equation on pressure";
        function squareReg = ThermoPower.Functions.squareReg;
      protected
        parameter Integer nXi = Medium.nXi "number of independent mass fractions";
        parameter Integer nX = Medium.nX "total number of mass fractions";
      public
        FlangeA infl(redeclare package Medium = Medium, m_flow(start = wnom, min = if allowFlowReversal then -Modelica.Constants.inf else 0));
        FlangeB outfl(redeclare package Medium = Medium, m_flow(start = -wnom, max = if allowFlowReversal then +Modelica.Constants.inf else 0));
      initial equation
        assert(wnom > 0, "Please set a positive value for wnom");
        assert(FFtype == FFtypes.NoFriction or dpnom > 0, "dpnom=0 not valid, it is also used in the homotopy trasformation during the inizialization");
        assert(not (FFtype == FFtypes.Kfnom and not Kfnom > 0), "Kfnom = 0 not valid, please set a positive value");
        assert(not (FFtype == FFtypes.OpPoint and not rhonom > 0), "rhonom = 0 not valid, please set a positive value");
        assert(not (FFtype == FFtypes.Cfnom and not Cfnom > 0), "Cfnom = 0 not valid, please set a positive value");
        assert(not (FFtype == FFtypes.Colebrook and not Dhyd > 0), "Dhyd = 0 not valid, please set a positive value");
        assert(not (FFtype == FFtypes.Colebrook and not e > 0), "e = 0 not valid, please set a positive value");
      end Flow1DBase;
    end BaseClasses;
  end Gas;

  package Thermal "Thermal models of heat transfer"
    extends Modelica.Icons.Package;

    connector DHTVolumes "Distributed Heat Terminal"
      parameter Integer N "Number of volumes";
      SI.Temperature[N] T "Temperature at the volumes";
      flow SI.Power[N] Q "Heat flow at the volumes";
    end DHTVolumes;

    model MetalTubeFV "Cylindrical metal tube model with Nw finite volumes"
      extends Icons.MetalWall;
      parameter Integer Nw = 1 "Number of volumes on the wall ports";
      parameter Integer Nt = 1 "Number of tubes in parallel";
      parameter SI.Length L "Tube length";
      parameter SI.Length rint "Internal radius (single tube)";
      parameter SI.Length rext "External radius (single tube)";
      parameter Real rhomcm "Metal heat capacity per unit volume [J/m^3.K]";
      parameter SI.ThermalConductivity lambda "Thermal conductivity";
      parameter Boolean WallRes = true "Wall thermal resistance accounted for";
      parameter SI.Temperature Tstartbar = 300 "Avarage temperature";
      parameter SI.Temperature Tstart1 = Tstartbar "Temperature start value - first volume";
      parameter SI.Temperature TstartN = Tstartbar "Temperature start value - last volume";
      parameter SI.Temperature[Nw] Tvolstart = Functions.linspaceExt(Tstart1, TstartN, Nw);
      parameter Choices.Init.Options initOpt = system.initOpt "Initialisation option";
      constant Real pi = Modelica.Constants.pi;
      final parameter SI.Area Am = (rext ^ 2 - rint ^ 2) * pi "Area of the metal tube cross-section";
      final parameter SI.HeatCapacity Cm = Nt * L * Am * rhomcm "Total heat capacity";
      outer ThermoPower.System system "System wide properties";
      SI.Temperature[Nw] Tvol(start = Tvolstart) "Volume temperatures";
      ThermoPower.Thermal.DHTVolumes int(final N = Nw, T(start = Tvolstart)) "Internal surface";
      ThermoPower.Thermal.DHTVolumes ext(final N = Nw, T(start = Tvolstart)) "External surface";
    initial equation
      if initOpt == Choices.Init.Options.noInit then
      elseif initOpt == Choices.Init.Options.fixedState then
        Tvol = Tvolstart;
      elseif initOpt == Choices.Init.Options.steadyState then
        der(Tvol) = zeros(Nw);
      elseif initOpt == Choices.Init.Options.steadyStateNoT then
      else
        assert(false, "Unsupported initialisation option");
      end if;
    equation
      assert(rext > rint, "External radius must be greater than internal radius");
      L / Nw * Nt * rhomcm * Am * der(Tvol) = int.Q + ext.Q "Energy balance";
      if WallRes then
        int.Q = lambda * (2 * pi * L / Nw) * (int.T - Tvol) / log((rint + rext) / (2 * rint)) * Nt "Heat conduction through the internal half-thickness";
        ext.Q = lambda * (2 * pi * L / Nw) * (ext.T - Tvol) / log(2 * rext / (rint + rext)) * Nt "Heat conduction through the external half-thickness";
      else
        ext.T = Tvol;
        int.T = Tvol;
      end if;
    end MetalTubeFV;

    model HeatExchangerTopologyFV "Connects two DHTVolumes ports according to a selected heat exchanger topology"
      extends Icons.HeatFlow;
      parameter Integer Nw "Number of volumes";
      replaceable model HeatExchangerTopology = HeatExchangerTopologies.CoCurrentFlow constrainedby ThermoPower.Thermal.BaseClasses.HeatExchangerTopologyData;
      HeatExchangerTopology HET(final Nw = Nw);
      Thermal.DHTVolumes side1(final N = Nw);
      Thermal.DHTVolumes side2(final N = Nw);
    equation
      for j in 1:Nw loop
        side2.T[HET.correspondingVolumes[j]] = side1.T[j];
        side2.Q[HET.correspondingVolumes[j]] + side1.Q[j] = 0;
      end for;
    end HeatExchangerTopologyFV;

    package HeatTransferFV "Heat transfer models for FV components"
      model IdealHeatTransfer "Delta T across the boundary layer is zero (infinite h.t.c.)"
        extends BaseClasses.DistributedHeatTransferFV(final useAverageTemperature = false);
      equation
        assert(Nw == Nf - 1, "Number of volumes Nw on wall side should be equal to number of volumes fluid side Nf - 1");
        for j in 1:Nw loop
          wall.T[j] = T[j + 1] "Ideal infinite heat transfer";
        end for;
      end IdealHeatTransfer;
    end HeatTransferFV;

    package HeatExchangerTopologies
      model CoCurrentFlow "Co-current flow"
        extends BaseClasses.HeatExchangerTopologyData(final correspondingVolumes = 1:Nw);
      end CoCurrentFlow;
    end HeatExchangerTopologies;

    package BaseClasses
      partial model DistributedHeatTransferFV "Base class for distributed heat transfer models - finite volumes"
        extends ThermoPower.Icons.HeatFlow;
        input Medium.ThermodynamicState[Nf] fluidState;
        input Medium.MassFlowRate[Nf] w;
        parameter Boolean useAverageTemperature = true "= true to use average temperature for heat transfer";
        ThermoPower.Thermal.DHTVolumes wall(final N = Nw);
        replaceable package Medium = Modelica.Media.Interfaces.PartialMedium "Medium model";
        parameter Integer Nf(min = 2) = 2 "Number of nodes on the fluid side";
        parameter Integer Nw = Nf - 1 "Number of volumes on the wall side";
        parameter Integer Nt(min = 1) "Number of tubes in parallel";
        parameter SI.Distance L "Tube length";
        parameter SI.Area A "Cross-sectional area (single tube)";
        parameter SI.Length omega "Wet perimeter of heat transfer surface (single tube)";
        parameter SI.Length Dhyd "Hydraulic Diameter (single tube)";
        parameter SI.MassFlowRate wnom "Nominal mass flow rate (single tube)";
        final parameter SI.Length l = L / Nw "Length of a single volume";
        parameter Modelica.SIunits.PerUnit kc = 1 "Correction factor for heat transfer";
        Medium.Temperature[Nf] T "Temperatures at the fluid side nodes";
        Medium.Temperature[Nw] Tw "Temperatures of the wall volumes";
        SI.Power[Nw] Qw "Heat flows entering from the wall volumes";
        SI.Power[Nf - 1] Qvol = Qw "Heat flows going to the fluid volumes";
        SI.Power Q "Total heat flow through lateral boundary";
      equation
        for j in 1:Nf loop
          T[j] = Medium.temperature(fluidState[j]);
        end for;
        Tw = wall.T;
        Qw = wall.Q;
        Q = sum(wall.Q);
      end DistributedHeatTransferFV;

      partial model HeatExchangerTopologyData "Base class for heat exchanger topology data"
        parameter Integer Nw "Number of volumes on both sides";
        parameter Integer[Nw] correspondingVolumes "Indeces of corresponding volumes";
      end HeatExchangerTopologyData;
    end BaseClasses;
  end Thermal;

  package Icons "Icons for ThermoPower library"
    extends Modelica.Icons.IconsPackage;

    partial model HeatFlow end HeatFlow;

    partial model MetalWall end MetalWall;

    package Gas "Icons for component using water/steam as working fluid"
      extends Modelica.Icons.Package;

      partial model SourceP end SourceP;

      partial model SourceW end SourceW;

      partial model Tube end Tube;

      partial model Compressor end Compressor;
    end Gas;
  end Icons;

  package Choices "Choice enumerations for ThermoPower models"
    extends Modelica.Icons.TypesPackage;

    package Flow1D
      type FFtypes = enumeration(Kfnom "Kfnom friction factor", OpPoint "Friction factor defined by operating point", Cfnom "Cfnom friction factor", Colebrook "Colebrook's equation", NoFriction "No friction") "Type, constants and menu choices to select the friction factor";
      type HCtypes = enumeration(Middle "Middle of the pipe", Upstream "At the inlet", Downstream "At the outlet") "Type, constants and menu choices to select the location of the hydraulic capacitance";
    end Flow1D;

    package TurboMachinery
      type TableTypes = enumeration(matrix "Explicitly supplied as parameter matrix table", file "Read from a file") "Type, constants and menu choices to select the representation of table matrix";
    end TurboMachinery;

    package Init "Options for initialisation"
      type Options = enumeration(noInit "No initial equations", fixedState "Fixed initial state variables", steadyState "Steady-state initialization", steadyStateNoP "Steady-state initialization except pressures (deprecated)", steadyStateNoT "Steady-state initialization except temperatures (deprecated)", steadyStateNoPT "Steady-state initialization except pressures and temperatures (deprecated)") "Type, constants and menu choices to select the initialisation options";
    end Init;
  end Choices;

  package Functions "Miscellaneous functions"
    extends Modelica.Icons.Package;

    function squareReg "Anti-symmetric square approximation with non-zero derivative in the origin"
      extends Modelica.Icons.Function;
      input Real x;
      input Real delta = 0.01 "Range of significant deviation from x^2*sgn(x)";
      output Real y;
    algorithm
      y := x * sqrt(x * x + delta * delta);
    end squareReg;

    function linspaceExt "Extended linspace function handling also the N=1 case"
      input Real x1;
      input Real x2;
      input Integer N;
      output Real[N] vec;
    algorithm
      vec := if N == 1 then {x1} else linspace(x1, x2, N);
    end linspaceExt;
  end Functions;

  package Units "Types with custom units"
    extends Modelica.Icons.Package;
    type HydraulicConductance = Real(final quantity = "HydraulicConductance", final unit = "(kg/s)/Pa");
    type HydraulicResistance = Real(final quantity = "HydraulicResistance", final unit = "Pa/(kg/s)");
  end Units;
  annotation(version = "3.1");
end ThermoPower;

package ModelicaServices "ModelicaServices (OpenModelica implementation) - Models and functions used in the Modelica Standard Library requiring a tool specific implementation"
  extends Modelica.Icons.Package;

  package Machine "Machine dependent constants"
    extends Modelica.Icons.Package;
    final constant Real eps = 1e-15 "Biggest number such that 1.0 + eps = 1.0";
    final constant Real small = 1e-60 "Smallest number such that small and -small are representable on the machine";
    final constant Real inf = 1e60 "Biggest Real number such that inf and -inf are representable on the machine";
    final constant Integer Integer_inf = OpenModelica.Internal.Architecture.integerMax() "Biggest Integer number such that Integer_inf and -Integer_inf are representable on the machine";
  end Machine;
  annotation(version = "4.0.0", versionDate = "2020-06-04", dateModified = "2020-06-04 11:00:00Z");
end ModelicaServices;

package Modelica "Modelica Standard Library - Version 3.2.3"
  extends Modelica.Icons.Package;

  package Blocks "Library of basic input/output control blocks (continuous, discrete, logical, table blocks)"
    import SI = Modelica.SIunits;
    extends Modelica.Icons.Package;

    package Interfaces "Library of connectors and partial models for input/output blocks"
      import Modelica.SIunits;
      extends Modelica.Icons.InterfacesPackage;
      connector RealInput = input Real "'input Real' as connector";
      connector RealOutput = output Real "'output Real' as connector";

      partial block SI2SO "2 Single Input / 1 Single Output continuous control block"
        extends Modelica.Blocks.Icons.Block;
        RealInput u1 "Connector of Real input signal 1";
        RealInput u2 "Connector of Real input signal 2";
        RealOutput y "Connector of Real output signal";
      end SI2SO;
    end Interfaces;

    package Tables "Library of blocks to interpolate in one and two-dimensional tables"
      extends Modelica.Icons.Package;

      block CombiTable2D "Table look-up in two dimensions (matrix/file)"
        extends Modelica.Blocks.Interfaces.SI2SO;
        extends Internal.CombiTable2DBase;
      equation
        if verboseExtrapolation and (extrapolation == Modelica.Blocks.Types.Extrapolation.LastTwoPoints or extrapolation == Modelica.Blocks.Types.Extrapolation.HoldLastPoint) then
          assert(noEvent(u1 >= u_min[1]), "
      Extrapolation warning: The value u1 (=" + String(u1) + ") must be greater or equal
      than the minimum abscissa value u_min[1] (=" + String(u_min[1]) + ") defined in the table.
          ", AssertionLevel.warning);
          assert(noEvent(u1 <= u_max[1]), "
      Extrapolation warning: The value u1 (=" + String(u1) + ") must be less or equal
      than the maximum abscissa value u_max[1] (=" + String(u_max[1]) + ") defined in the table.
          ", AssertionLevel.warning);
          assert(noEvent(u2 >= u_min[2]), "
      Extrapolation warning: The value u2 (=" + String(u2) + ") must be greater or equal
      than the minimum abscissa value u_min[2] (=" + String(u_min[2]) + ") defined in the table.
          ", AssertionLevel.warning);
          assert(noEvent(u2 <= u_max[2]), "
      Extrapolation warning: The value u2 (=" + String(u2) + ") must be less or equal
      than the maximum abscissa value u_max[2] (=" + String(u_max[2]) + ") defined in the table.
          ", AssertionLevel.warning);
        end if;
        if smoothness == Modelica.Blocks.Types.Smoothness.ConstantSegments then
          y = Internal.getTable2DValueNoDer(tableID, u1, u2);
        else
          y = Internal.getTable2DValue(tableID, u1, u2);
        end if;
      end CombiTable2D;

      package Internal "Internal external object definitions for table functions that should not be directly utilized by the user"
        extends Modelica.Icons.InternalPackage;

        partial block CombiTable2DBase "Base class for variants of CombiTable2D"
          parameter Boolean tableOnFile = false "= true, if table is defined on file or in function usertab";
          parameter Real[:, :] table = fill(0.0, 0, 2) "Table matrix (grid u1 = first column, grid u2 = first row; e.g., table=[0, 0; 0, 1])";
          parameter String tableName = "NoName" "Table name on file or in function usertab (see docu)";
          parameter String fileName = "NoName" "File where matrix is stored";
          parameter Boolean verboseRead = true "= true, if info message that file is loading is to be printed";
          parameter Modelica.Blocks.Types.Smoothness smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments "Smoothness of table interpolation";
          parameter Modelica.Blocks.Types.Extrapolation extrapolation = Modelica.Blocks.Types.Extrapolation.LastTwoPoints "Extrapolation of data outside the definition range";
          parameter Boolean verboseExtrapolation = false "= true, if warning messages are to be printed if table input is outside the definition range";
          final parameter Real[2] u_min = getTable2DAbscissaUmin(tableID) "Minimum abscissa value defined in table";
          final parameter Real[2] u_max = getTable2DAbscissaUmax(tableID) "Maximum abscissa value defined in table";
        protected
          parameter Modelica.Blocks.Types.ExternalCombiTable2D tableID = Modelica.Blocks.Types.ExternalCombiTable2D(if tableOnFile then tableName else "NoName", if tableOnFile and fileName <> "NoName" and not Modelica.Utilities.Strings.isEmpty(fileName) then fileName else "NoName", table, smoothness, extrapolation, if tableOnFile then verboseRead else false) "External table object";
        equation
          if tableOnFile then
            assert(tableName <> "NoName", "tableOnFile = true and no table name given");
          else
            assert(size(table, 1) > 0 and size(table, 2) > 0, "tableOnFile = false and parameter table is an empty matrix");
          end if;
        end CombiTable2DBase;

        function getTable2DValue "Interpolate 2-dim. table defined by matrix"
          extends Modelica.Icons.Function;
          input Modelica.Blocks.Types.ExternalCombiTable2D tableID;
          input Real u1;
          input Real u2;
          output Real y;
          external "C" y = ModelicaStandardTables_CombiTable2D_getValue(tableID, u1, u2) annotation(Library = {"ModelicaStandardTables", "ModelicaIO", "ModelicaMatIO", "zlib"}, derivative = getDerTable2DValue);
          annotation(derivative = getDerTable2DValue);
        end getTable2DValue;

        function getTable2DValueNoDer "Interpolate 2-dim. table defined by matrix (but do not provide a derivative function)"
          extends Modelica.Icons.Function;
          input Modelica.Blocks.Types.ExternalCombiTable2D tableID;
          input Real u1;
          input Real u2;
          output Real y;
          external "C" y = ModelicaStandardTables_CombiTable2D_getValue(tableID, u1, u2) annotation(Library = {"ModelicaStandardTables", "ModelicaIO", "ModelicaMatIO", "zlib"});
        end getTable2DValueNoDer;

        function getDerTable2DValue "Derivative of interpolated 2-dim. table defined by matrix"
          extends Modelica.Icons.Function;
          input Modelica.Blocks.Types.ExternalCombiTable2D tableID;
          input Real u1;
          input Real u2;
          input Real der_u1;
          input Real der_u2;
          output Real der_y;
          external "C" der_y = ModelicaStandardTables_CombiTable2D_getDerValue(tableID, u1, u2, der_u1, der_u2) annotation(Library = {"ModelicaStandardTables", "ModelicaIO", "ModelicaMatIO", "zlib"});
        end getDerTable2DValue;

        function getTable2DAbscissaUmin "Return minimum abscissa value of 2-dim. table defined by matrix"
          extends Modelica.Icons.Function;
          input Modelica.Blocks.Types.ExternalCombiTable2D tableID;
          output Real[2] uMin "Minimum abscissa value in table";
          external "C" ModelicaStandardTables_CombiTable2D_minimumAbscissa(tableID, uMin) annotation(Library = {"ModelicaStandardTables", "ModelicaIO", "ModelicaMatIO", "zlib"});
        end getTable2DAbscissaUmin;

        function getTable2DAbscissaUmax "Return maximum abscissa value of 2-dim. table defined by matrix"
          extends Modelica.Icons.Function;
          input Modelica.Blocks.Types.ExternalCombiTable2D tableID;
          output Real[2] uMax "Maximum abscissa value in table";
          external "C" ModelicaStandardTables_CombiTable2D_maximumAbscissa(tableID, uMax) annotation(Library = {"ModelicaStandardTables", "ModelicaIO", "ModelicaMatIO", "zlib"});
        end getTable2DAbscissaUmax;
      end Internal;
    end Tables;

    package Types "Library of constants, external objects and types with choices, especially to build menus"
      extends Modelica.Icons.TypesPackage;
      type Smoothness = enumeration(LinearSegments "Table points are linearly interpolated", ContinuousDerivative "Table points are interpolated (by Akima splines) such that the first derivative is continuous", ConstantSegments "Table points are not interpolated, but the value from the previous abscissa point is returned", MonotoneContinuousDerivative1 "Table points are interpolated (by Fritsch-Butland splines) such that the monotonicity is preserved and the first derivative is continuous", MonotoneContinuousDerivative2 "Table points are interpolated (by Steffen splines) such that the monotonicity is preserved and the first derivative is continuous") "Enumeration defining the smoothness of table interpolation";
      type Extrapolation = enumeration(HoldLastPoint "Hold the first/last table point outside of the table scope", LastTwoPoints "Extrapolate by using the derivative at the first/last table points outside of the table scope", Periodic "Repeat the table scope periodically", NoExtrapolation "Extrapolation triggers an error") "Enumeration defining the extrapolation of table interpolation";

      class ExternalCombiTable2D "External object of 2-dim. table defined by matrix"
        extends ExternalObject;

        function constructor "Initialize 2-dim. table defined by matrix"
          extends Modelica.Icons.Function;
          input String tableName "Table name";
          input String fileName "File name";
          input Real[:, :] table;
          input Modelica.Blocks.Types.Smoothness smoothness;
          input Modelica.Blocks.Types.Extrapolation extrapolation = Modelica.Blocks.Types.Extrapolation.LastTwoPoints;
          input Boolean verboseRead = true "= true: Print info message; = false: No info message";
          output ExternalCombiTable2D externalCombiTable2D;
          external "C" externalCombiTable2D = ModelicaStandardTables_CombiTable2D_init2(fileName, tableName, table, size(table, 1), size(table, 2), smoothness, extrapolation, verboseRead) annotation(Library = {"ModelicaStandardTables", "ModelicaIO", "ModelicaMatIO", "zlib"});
        end constructor;

        function destructor "Terminate 2-dim. table defined by matrix"
          extends Modelica.Icons.Function;
          input ExternalCombiTable2D externalCombiTable2D;
          external "C" ModelicaStandardTables_CombiTable2D_close(externalCombiTable2D) annotation(Library = {"ModelicaStandardTables", "ModelicaIO", "ModelicaMatIO", "zlib"});
        end destructor;
      end ExternalCombiTable2D;
    end Types;

    package Icons "Icons for Blocks"
      extends Modelica.Icons.IconsPackage;

      partial block Block "Basic graphical layout of input/output block" end Block;
    end Icons;
  end Blocks;

  package Mechanics "Library of 1-dim. and 3-dim. mechanical components (multi-body, rotational, translational)"
    extends Modelica.Icons.Package;

    package Rotational "Library to model 1-dimensional, rotational mechanical systems"
      extends Modelica.Icons.Package;
      import SI = Modelica.SIunits;

      package Interfaces "Connectors and partial models for 1D rotational mechanical components"
        extends Modelica.Icons.InterfacesPackage;

        connector Flange "One-dimensional rotational flange"
          SI.Angle phi "Absolute rotation angle of flange";
          flow SI.Torque tau "Cut torque in the flange";
        end Flange;

        connector Flange_a "One-dimensional rotational flange of a shaft (filled circle icon)"
          extends Flange;
        end Flange_a;

        connector Flange_b "One-dimensional rotational flange of a shaft (non-filled circle icon)"
          extends Flange;
        end Flange_b;
      end Interfaces;
    end Rotational;
  end Mechanics;

  package Media "Library of media property models"
    extends Modelica.Icons.Package;
    import SI = Modelica.SIunits;
    import Cv = Modelica.SIunits.Conversions;

    package Interfaces "Interfaces for media models"
      extends Modelica.Icons.InterfacesPackage;

      partial package PartialMedium "Partial medium properties (base package of all media packages)"
        extends Modelica.Media.Interfaces.Types;
        extends Modelica.Icons.MaterialPropertiesPackage;
        constant Modelica.Media.Interfaces.Choices.IndependentVariables ThermoStates "Enumeration type for independent variables";
        constant String mediumName = "unusablePartialMedium" "Name of the medium";
        constant String[:] substanceNames = {mediumName} "Names of the mixture substances. Set substanceNames={mediumName} if only one substance.";
        constant String[:] extraPropertiesNames = fill("", 0) "Names of the additional (extra) transported properties. Set extraPropertiesNames=fill(\"\",0) if unused";
        constant Boolean singleState "= true, if u and d are not a function of pressure";
        constant Boolean reducedX = true "= true if medium contains the equation sum(X) = 1.0; set reducedX=true if only one substance (see docu for details)";
        constant Boolean fixedX = false "= true if medium contains the equation X = reference_X";
        constant AbsolutePressure reference_p = 101325 "Reference pressure of Medium: default 1 atmosphere";
        constant MassFraction[nX] reference_X = fill(1 / nX, nX) "Default mass fractions of medium";
        constant AbsolutePressure p_default = 101325 "Default value for pressure of medium (for initialization)";
        constant Temperature T_default = Modelica.SIunits.Conversions.from_degC(20) "Default value for temperature of medium (for initialization)";
        constant MassFraction[nX] X_default = reference_X "Default value for mass fractions of medium (for initialization)";
        final constant Integer nS = size(substanceNames, 1) "Number of substances";
        constant Integer nX = nS "Number of mass fractions";
        constant Integer nXi = if fixedX then 0 else if reducedX then nS - 1 else nS "Number of structurally independent mass fractions (see docu for details)";
        final constant Integer nC = size(extraPropertiesNames, 1) "Number of extra (outside of standard mass-balance) transported properties";
        replaceable record FluidConstants = Modelica.Media.Interfaces.Types.Basic.FluidConstants "Critical, triple, molecular and other standard data of fluid";

        replaceable record ThermodynamicState "Minimal variable set that is available as input argument to every medium function"
          extends Modelica.Icons.Record;
        end ThermodynamicState;

        replaceable partial model BaseProperties "Base properties (p, d, T, h, u, R, MM and, if applicable, X and Xi) of a medium"
          InputAbsolutePressure p "Absolute pressure of medium";
          InputMassFraction[nXi] Xi(start = reference_X[1:nXi]) "Structurally independent mass fractions";
          InputSpecificEnthalpy h "Specific enthalpy of medium";
          Density d "Density of medium";
          Temperature T "Temperature of medium";
          MassFraction[nX] X(start = reference_X) "Mass fractions (= (component mass)/total mass  m_i/m)";
          SpecificInternalEnergy u "Specific internal energy of medium";
          SpecificHeatCapacity R "Gas constant (of mixture if applicable)";
          MolarMass MM "Molar mass (of mixture or single fluid)";
          ThermodynamicState state "Thermodynamic state record for optional functions";
          parameter Boolean preferredMediumStates = false "= true if StateSelect.prefer shall be used for the independent property variables of the medium" annotation(Evaluate = true);
          parameter Boolean standardOrderComponents = true "If true, and reducedX = true, the last element of X will be computed from the other ones";
          SI.Conversions.NonSIunits.Temperature_degC T_degC = Modelica.SIunits.Conversions.to_degC(T) "Temperature of medium in [degC]";
          SI.Conversions.NonSIunits.Pressure_bar p_bar = Modelica.SIunits.Conversions.to_bar(p) "Absolute pressure of medium in [bar]";
          connector InputAbsolutePressure = input SI.AbsolutePressure "Pressure as input signal connector";
          connector InputSpecificEnthalpy = input SI.SpecificEnthalpy "Specific enthalpy as input signal connector";
          connector InputMassFraction = input SI.MassFraction "Mass fraction as input signal connector";
        equation
          if standardOrderComponents then
            Xi = X[1:nXi];
            if fixedX then
              X = reference_X;
            end if;
            if reducedX and not fixedX then
              X[nX] = 1 - sum(Xi);
            end if;
            for i in 1:nX loop
              assert(X[i] >= (-1.e-5) and X[i] <= 1 + 1.e-5, "Mass fraction X[" + String(i) + "] = " + String(X[i]) + "of substance " + substanceNames[i] + "\nof medium " + mediumName + " is not in the range 0..1");
            end for;
          end if;
          assert(p >= 0.0, "Pressure (= " + String(p) + " Pa) of medium \"" + mediumName + "\" is negative\n(Temperature = " + String(T) + " K)");
        end BaseProperties;

        replaceable partial function setState_pTX "Return thermodynamic state as function of p, T and composition X or Xi"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input Temperature T "Temperature";
          input MassFraction[:] X = reference_X "Mass fractions";
          output ThermodynamicState state "Thermodynamic state record";
        end setState_pTX;

        replaceable partial function setState_phX "Return thermodynamic state as function of p, h and composition X or Xi"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SpecificEnthalpy h "Specific enthalpy";
          input MassFraction[:] X = reference_X "Mass fractions";
          output ThermodynamicState state "Thermodynamic state record";
        end setState_phX;

        replaceable partial function setState_psX "Return thermodynamic state as function of p, s and composition X or Xi"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SpecificEntropy s "Specific entropy";
          input MassFraction[:] X = reference_X "Mass fractions";
          output ThermodynamicState state "Thermodynamic state record";
        end setState_psX;

        replaceable partial function setState_dTX "Return thermodynamic state as function of d, T and composition X or Xi"
          extends Modelica.Icons.Function;
          input Density d "Density";
          input Temperature T "Temperature";
          input MassFraction[:] X = reference_X "Mass fractions";
          output ThermodynamicState state "Thermodynamic state record";
        end setState_dTX;

        replaceable partial function setSmoothState "Return thermodynamic state so that it smoothly approximates: if x > 0 then state_a else state_b"
          extends Modelica.Icons.Function;
          input Real x "m_flow or dp";
          input ThermodynamicState state_a "Thermodynamic state if x > 0";
          input ThermodynamicState state_b "Thermodynamic state if x < 0";
          input Real x_small(min = 0) "Smooth transition in the region -x_small < x < x_small";
          output ThermodynamicState state "Smooth thermodynamic state for all x (continuous and differentiable)";
        end setSmoothState;

        replaceable partial function dynamicViscosity "Return dynamic viscosity"
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output DynamicViscosity eta "Dynamic viscosity";
        end dynamicViscosity;

        replaceable partial function thermalConductivity "Return thermal conductivity"
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output ThermalConductivity lambda "Thermal conductivity";
        end thermalConductivity;

        replaceable partial function pressure "Return pressure"
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output AbsolutePressure p "Pressure";
        end pressure;

        replaceable partial function temperature "Return temperature"
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output Temperature T "Temperature";
        end temperature;

        replaceable partial function density "Return density"
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output Density d "Density";
        end density;

        replaceable partial function specificEnthalpy "Return specific enthalpy"
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output SpecificEnthalpy h "Specific enthalpy";
        end specificEnthalpy;

        replaceable partial function specificInternalEnergy "Return specific internal energy"
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output SpecificEnergy u "Specific internal energy";
        end specificInternalEnergy;

        replaceable partial function specificEntropy "Return specific entropy"
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output SpecificEntropy s "Specific entropy";
        end specificEntropy;

        replaceable partial function specificGibbsEnergy "Return specific Gibbs energy"
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output SpecificEnergy g "Specific Gibbs energy";
        end specificGibbsEnergy;

        replaceable partial function specificHelmholtzEnergy "Return specific Helmholtz energy"
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output SpecificEnergy f "Specific Helmholtz energy";
        end specificHelmholtzEnergy;

        replaceable partial function specificHeatCapacityCp "Return specific heat capacity at constant pressure"
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output SpecificHeatCapacity cp "Specific heat capacity at constant pressure";
        end specificHeatCapacityCp;

        replaceable partial function specificHeatCapacityCv "Return specific heat capacity at constant volume"
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output SpecificHeatCapacity cv "Specific heat capacity at constant volume";
        end specificHeatCapacityCv;

        function heatCapacity_cv = specificHeatCapacityCv "Alias for deprecated name";

        replaceable partial function isentropicExponent "Return isentropic exponent"
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output IsentropicExponent gamma "Isentropic exponent";
        end isentropicExponent;

        replaceable partial function isentropicEnthalpy "Return isentropic enthalpy"
          extends Modelica.Icons.Function;
          input AbsolutePressure p_downstream "Downstream pressure";
          input ThermodynamicState refState "Reference state for entropy";
          output SpecificEnthalpy h_is "Isentropic enthalpy";
        end isentropicEnthalpy;

        replaceable partial function velocityOfSound "Return velocity of sound"
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output VelocityOfSound a "Velocity of sound";
        end velocityOfSound;

        replaceable partial function isobaricExpansionCoefficient "Return overall the isobaric expansion coefficient beta"
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output IsobaricExpansionCoefficient beta "Isobaric expansion coefficient";
        end isobaricExpansionCoefficient;

        replaceable partial function isothermalCompressibility "Return overall the isothermal compressibility factor"
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output SI.IsothermalCompressibility kappa "Isothermal compressibility";
        end isothermalCompressibility;

        replaceable partial function density_derp_T "Return density derivative w.r.t. pressure at const temperature"
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output DerDensityByPressure ddpT "Density derivative w.r.t. pressure";
        end density_derp_T;

        replaceable partial function density_derT_p "Return density derivative w.r.t. temperature at constant pressure"
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output DerDensityByTemperature ddTp "Density derivative w.r.t. temperature";
        end density_derT_p;

        replaceable partial function density_derX "Return density derivative w.r.t. mass fraction"
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output Density[nX] dddX "Derivative of density w.r.t. mass fraction";
        end density_derX;

        replaceable partial function molarMass "Return the molar mass of the medium"
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output MolarMass MM "Mixture molar mass";
        end molarMass;

        replaceable function specificEnthalpy_pTX "Return specific enthalpy from p, T, and X or Xi"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input Temperature T "Temperature";
          input MassFraction[:] X = reference_X "Mass fractions";
          output SpecificEnthalpy h "Specific enthalpy";
        algorithm
          h := specificEnthalpy(setState_pTX(p, T, X));
          annotation(inverse(T = temperature_phX(p, h, X)));
        end specificEnthalpy_pTX;

        replaceable function temperature_phX "Return temperature from p, h, and X or Xi"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SpecificEnthalpy h "Specific enthalpy";
          input MassFraction[:] X = reference_X "Mass fractions";
          output Temperature T "Temperature";
        algorithm
          T := temperature(setState_phX(p, h, X));
        end temperature_phX;

        type MassFlowRate = SI.MassFlowRate(quantity = "MassFlowRate." + mediumName, min = -1.0e5, max = 1.e5) "Type for mass flow rate with medium specific attributes";
      end PartialMedium;

      partial package PartialPureSubstance "Base class for pure substances of one chemical substance"
        extends PartialMedium(final reducedX = true, final fixedX = true);

        redeclare replaceable partial model extends BaseProperties(final standardOrderComponents = true) end BaseProperties;
      end PartialPureSubstance;

      package Choices "Types, constants to define menu choices"
        extends Modelica.Icons.Package;
        type IndependentVariables = enumeration(T "Temperature", pT "Pressure, Temperature", ph "Pressure, Specific Enthalpy", phX "Pressure, Specific Enthalpy, Mass Fraction", pTX "Pressure, Temperature, Mass Fractions", dTX "Density, Temperature, Mass Fractions") "Enumeration defining the independent variables of a medium";
        type ReferenceEnthalpy = enumeration(ZeroAt0K "The enthalpy is 0 at 0 K (default), if the enthalpy of formation is excluded", ZeroAt25C "The enthalpy is 0 at 25 degC, if the enthalpy of formation is excluded", UserDefined "The user-defined reference enthalpy is used at 293.15 K (25 degC)") "Enumeration defining the reference enthalpy of a medium" annotation(Evaluate = true);
      end Choices;

      package Types "Types to be used in fluid models"
        extends Modelica.Icons.Package;
        type AbsolutePressure = SI.AbsolutePressure(min = 0, max = 1.e8, nominal = 1.e5, start = 1.e5) "Type for absolute pressure with medium specific attributes";
        type Density = SI.Density(min = 0, max = 1.e5, nominal = 1, start = 1) "Type for density with medium specific attributes";
        type DynamicViscosity = SI.DynamicViscosity(min = 0, max = 1.e8, nominal = 1.e-3, start = 1.e-3) "Type for dynamic viscosity with medium specific attributes";
        type MassFraction = Real(quantity = "MassFraction", final unit = "kg/kg", min = 0, max = 1, nominal = 0.1) "Type for mass fraction with medium specific attributes";
        type MolarMass = SI.MolarMass(min = 0.001, max = 0.25, nominal = 0.032) "Type for molar mass with medium specific attributes";
        type MolarVolume = SI.MolarVolume(min = 1e-6, max = 1.0e6, nominal = 1.0) "Type for molar volume with medium specific attributes";
        type IsentropicExponent = SI.RatioOfSpecificHeatCapacities(min = 1, max = 500000, nominal = 1.2, start = 1.2) "Type for isentropic exponent with medium specific attributes";
        type SpecificEnergy = SI.SpecificEnergy(min = -1.0e8, max = 1.e8, nominal = 1.e6) "Type for specific energy with medium specific attributes";
        type SpecificInternalEnergy = SpecificEnergy "Type for specific internal energy with medium specific attributes";
        type SpecificEnthalpy = SI.SpecificEnthalpy(min = -1.0e10, max = 1.e10, nominal = 1.e6) "Type for specific enthalpy with medium specific attributes";
        type SpecificEntropy = SI.SpecificEntropy(min = -1.e7, max = 1.e7, nominal = 1.e3) "Type for specific entropy with medium specific attributes";
        type SpecificHeatCapacity = SI.SpecificHeatCapacity(min = 0, max = 1.e7, nominal = 1.e3, start = 1.e3) "Type for specific heat capacity with medium specific attributes";
        type Temperature = SI.Temperature(min = 1, max = 1.e4, nominal = 300, start = 288.15) "Type for temperature with medium specific attributes";
        type ThermalConductivity = SI.ThermalConductivity(min = 0, max = 500, nominal = 1, start = 1) "Type for thermal conductivity with medium specific attributes";
        type VelocityOfSound = SI.Velocity(min = 0, max = 1.e5, nominal = 1000, start = 1000) "Type for velocity of sound with medium specific attributes";
        type ExtraProperty = Real(min = 0.0, start = 1.0) "Type for unspecified, mass-specific property transported by flow";
        type IsobaricExpansionCoefficient = Real(min = 0, max = 1.0e8, unit = "1/K") "Type for isobaric expansion coefficient with medium specific attributes";
        type DipoleMoment = Real(min = 0.0, max = 2.0, unit = "debye", quantity = "ElectricDipoleMoment") "Type for dipole moment with medium specific attributes";
        type DerDensityByPressure = SI.DerDensityByPressure "Type for partial derivative of density with respect to pressure with medium specific attributes";
        type DerDensityByTemperature = SI.DerDensityByTemperature "Type for partial derivative of density with respect to temperature with medium specific attributes";

        package Basic "The most basic version of a record used in several degrees of detail"
          extends Icons.Package;

          record FluidConstants "Critical, triple, molecular and other standard data of fluid"
            extends Modelica.Icons.Record;
            String iupacName "Complete IUPAC name (or common name, if non-existent)";
            String casRegistryNumber "Chemical abstracts sequencing number (if it exists)";
            String chemicalFormula "Chemical formula, (brutto, nomenclature according to Hill";
            String structureFormula "Chemical structure formula";
            MolarMass molarMass "Molar mass";
          end FluidConstants;
        end Basic;

        package IdealGas "The ideal gas version of a record used in several degrees of detail"
          extends Icons.Package;

          record FluidConstants "Extended fluid constants"
            extends Modelica.Media.Interfaces.Types.Basic.FluidConstants;
            Temperature criticalTemperature "Critical temperature";
            AbsolutePressure criticalPressure "Critical pressure";
            MolarVolume criticalMolarVolume "Critical molar Volume";
            Real acentricFactor "Pitzer acentric factor";
            Temperature meltingPoint "Melting point at 101325 Pa";
            Temperature normalBoilingPoint "Normal boiling point (at 101325 Pa)";
            DipoleMoment dipoleMoment "Dipole moment of molecule in Debye (1 debye = 3.33564e10-30 C.m)";
            Boolean hasIdealGasHeatCapacity = false "True if ideal gas heat capacity is available";
            Boolean hasCriticalData = false "True if critical data are known";
            Boolean hasDipoleMoment = false "True if a dipole moment known";
            Boolean hasFundamentalEquation = false "True if a fundamental equation";
            Boolean hasLiquidHeatCapacity = false "True if liquid heat capacity is available";
            Boolean hasSolidHeatCapacity = false "True if solid heat capacity is available";
            Boolean hasAccurateViscosityData = false "True if accurate data for a viscosity function is available";
            Boolean hasAccurateConductivityData = false "True if accurate data for thermal conductivity is available";
            Boolean hasVapourPressureCurve = false "True if vapour pressure data, e.g., Antoine coefficients are known";
            Boolean hasAcentricFactor = false "True if Pitzer acentric factor is known";
            SpecificEnthalpy HCRIT0 = 0.0 "Critical specific enthalpy of the fundamental equation";
            SpecificEntropy SCRIT0 = 0.0 "Critical specific entropy of the fundamental equation";
            SpecificEnthalpy deltah = 0.0 "Difference between specific enthalpy model (h_m) and f.eq. (h_f) (h_m - h_f)";
            SpecificEntropy deltas = 0.0 "Difference between specific enthalpy model (s_m) and f.eq. (s_f) (s_m - s_f)";
          end FluidConstants;
        end IdealGas;
      end Types;
    end Interfaces;

    package Common "Data structures and fundamental functions for fluid properties"
      extends Modelica.Icons.Package;
      constant Real MINPOS = 1.0e-9 "Minimal value for physical variables which are always > 0.0";

      function smoothStep "Approximation of a general step, such that the characteristic is continuous and differentiable"
        extends Modelica.Icons.Function;
        input Real x "Abscissa value";
        input Real y1 "Ordinate value for x > 0";
        input Real y2 "Ordinate value for x < 0";
        input Real x_small(min = 0) = 1e-5 "Approximation of step for -x_small <= x <= x_small; x_small > 0 required";
        output Real y "Ordinate value to approximate y = if x > 0 then y1 else y2";
      algorithm
        y := smooth(1, if x > x_small then y1 else if x < (-x_small) then y2 else if abs(x_small) > 0 then x / x_small * ((x / x_small) ^ 2 - 3) * (y2 - y1) / 4 + (y1 + y2) / 2 else (y1 + y2) / 2);
        annotation(Inline = true, smoothOrder = 1);
      end smoothStep;

      package OneNonLinearEquation "Determine solution of a non-linear algebraic equation in one unknown without derivatives in a reliable and efficient way"
        extends Modelica.Icons.Package;

        replaceable record f_nonlinear_Data "Data specific for function f_nonlinear"
          extends Modelica.Icons.Record;
        end f_nonlinear_Data;

        replaceable partial function f_nonlinear "Nonlinear algebraic equation in one unknown: y = f_nonlinear(x,p,X)"
          extends Modelica.Icons.Function;
          input Real x "Independent variable of function";
          input Real p = 0.0 "Disregarded variables (here always used for pressure)";
          input Real[:] X = fill(0, 0) "Disregarded variables (her always used for composition)";
          input f_nonlinear_Data f_nonlinear_data "Additional data for the function";
          output Real y "= f_nonlinear(x)";
        end f_nonlinear;

        replaceable function solve "Solve f_nonlinear(x_zero)=y_zero; f_nonlinear(x_min) - y_zero and f_nonlinear(x_max)-y_zero must have different sign"
          import Modelica.Utilities.Streams.error;
          extends Modelica.Icons.Function;
          input Real y_zero "Determine x_zero, such that f_nonlinear(x_zero) = y_zero";
          input Real x_min "Minimum value of x";
          input Real x_max "Maximum value of x";
          input Real pressure = 0.0 "Disregarded variables (here always used for pressure)";
          input Real[:] X = fill(0, 0) "Disregarded variables (here always used for composition)";
          input f_nonlinear_Data f_nonlinear_data "Additional data for function f_nonlinear";
          input Real x_tol = 100 * Modelica.Constants.eps "Relative tolerance of the result";
          output Real x_zero "f_nonlinear(x_zero) = y_zero";
        protected
          constant Real eps = Modelica.Constants.eps "Machine epsilon";
          constant Real x_eps = 1e-10 "Slight modification of x_min, x_max, since x_min, x_max are usually exactly at the borders T_min/h_min and then small numeric noise may make the interval invalid";
          Real x_min2 = x_min - x_eps;
          Real x_max2 = x_max + x_eps;
          Real a = x_min2 "Current best minimum interval value";
          Real b = x_max2 "Current best maximum interval value";
          Real c "Intermediate point a <= c <= b";
          Real d;
          Real e "b - a";
          Real m;
          Real s;
          Real p;
          Real q;
          Real r;
          Real tol;
          Real fa "= f_nonlinear(a) - y_zero";
          Real fb "= f_nonlinear(b) - y_zero";
          Real fc;
          Boolean found = false;
        algorithm
          fa := f_nonlinear(x_min2, pressure, X, f_nonlinear_data) - y_zero;
          fb := f_nonlinear(x_max2, pressure, X, f_nonlinear_data) - y_zero;
          fc := fb;
          if fa > 0.0 and fb > 0.0 or fa < 0.0 and fb < 0.0 then
            error("The arguments x_min and x_max to OneNonLinearEquation.solve(..)\n" + "do not bracket the root of the single non-linear equation:\n" + "  x_min  = " + String(x_min2) + "\n" + "  x_max  = " + String(x_max2) + "\n" + "  y_zero = " + String(y_zero) + "\n" + "  fa = f(x_min) - y_zero = " + String(fa) + "\n" + "  fb = f(x_max) - y_zero = " + String(fb) + "\n" + "fa and fb must have opposite sign which is not the case");
          else
          end if;
          c := a;
          fc := fa;
          e := b - a;
          d := e;
          while not found loop
            if abs(fc) < abs(fb) then
              a := b;
              b := c;
              c := a;
              fa := fb;
              fb := fc;
              fc := fa;
            else
            end if;
            tol := 2 * eps * abs(b) + x_tol;
            m := (c - b) / 2;
            if abs(m) <= tol or fb == 0.0 then
              found := true;
              x_zero := b;
            else
              if abs(e) < tol or abs(fa) <= abs(fb) then
                e := m;
                d := e;
              else
                s := fb / fa;
                if a == c then
                  p := 2 * m * s;
                  q := 1 - s;
                else
                  q := fa / fc;
                  r := fb / fc;
                  p := s * (2 * m * q * (q - r) - (b - a) * (r - 1));
                  q := (q - 1) * (r - 1) * (s - 1);
                end if;
                if p > 0 then
                  q := -q;
                else
                  p := -p;
                end if;
                s := e;
                e := d;
                if 2 * p < 3 * m * q - abs(tol * q) and p < abs(0.5 * s * q) then
                  d := p / q;
                else
                  e := m;
                  d := e;
                end if;
              end if;
              a := b;
              fa := fb;
              b := b + (if abs(d) > tol then d else if m > 0 then tol else -tol);
              fb := f_nonlinear(b, pressure, X, f_nonlinear_data) - y_zero;
              if fb > 0 and fc > 0 or fb < 0 and fc < 0 then
                c := a;
                fc := fa;
                e := b - a;
                d := e;
              else
              end if;
            end if;
          end while;
        end solve;
      end OneNonLinearEquation;
    end Common;

    package Air "Medium models for air"
      extends Modelica.Icons.VariantsPackage;

      package DryAirNasa "Air: Detailed dry air model as ideal gas (200..6000 K)"
        extends Modelica.Icons.MaterialProperty;
        extends IdealGases.Common.SingleGasNasa(mediumName = "Air", data = IdealGases.Common.SingleGasesData.Air, fluidConstants = {IdealGases.Common.FluidData.N2});

        redeclare function dynamicViscosity "Return dynamic viscosity of dry air (simple polynomial, moisture influence small, valid from 123.15 K to 1273.15 K, outside of this range linear extrapolation is used)"
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output DynamicViscosity eta "Dynamic viscosity";
          import Modelica.Media.Incompressible.TableBased.Polynomials_Temp;
        algorithm
          eta := 1e-6 * Polynomials_Temp.evaluateWithRange({9.7391102886305869E-15, -3.1353724870333906E-11, 4.3004876595642225E-08, -3.8228016291758240E-05, 5.0427874367180762E-02, 1.7239260139242528E+01}, Cv.to_degC(123.15), Cv.to_degC(1273.15), Cv.to_degC(state.T));
          annotation(smoothOrder = 2);
        end dynamicViscosity;

        redeclare function thermalConductivity "Return thermal conductivity of dry air (simple polynomial, moisture influence small, valid from 123.15 K to 1273.15 K, outside of this range linear extrapolation is used)"
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          input Integer method = 1 "Dummy for compatibility reasons";
          output ThermalConductivity lambda "Thermal conductivity";
          import Modelica.Media.Incompressible.TableBased.Polynomials_Temp;
          import Cv = Modelica.SIunits.Conversions;
        algorithm
          lambda := 1e-3 * Polynomials_Temp.evaluateWithRange({6.5691470817717812E-15, -3.4025961923050509E-11, 5.3279284846303157E-08, -4.5340839289219472E-05, 7.6129675309037664E-02, 2.4169481088097051E+01}, Cv.to_degC(123.15), Cv.to_degC(1273.15), Cv.to_degC(state.T));
          annotation(smoothOrder = 2);
        end thermalConductivity;
      end DryAirNasa;
    end Air;

    package IdealGases "Data and models of ideal gases (single, fixed and dynamic mixtures) from NASA source"
      extends Modelica.Icons.VariantsPackage;

      package Common "Common packages and data for the ideal gas models"
        extends Modelica.Icons.Package;

        record DataRecord "Coefficient data record for properties of ideal gases based on NASA source"
          extends Modelica.Icons.Record;
          String name "Name of ideal gas";
          SI.MolarMass MM "Molar mass";
          SI.SpecificEnthalpy Hf "Enthalpy of formation at 298.15K";
          SI.SpecificEnthalpy H0 "H0(298.15K) - H0(0K)";
          SI.Temperature Tlimit "Temperature limit between low and high data sets";
          Real[7] alow "Low temperature coefficients a";
          Real[2] blow "Low temperature constants b";
          Real[7] ahigh "High temperature coefficients a";
          Real[2] bhigh "High temperature constants b";
          SI.SpecificHeatCapacity R "Gas constant";
        end DataRecord;

        package Functions "Basic Functions for ideal gases: cp, h, s, thermal conductivity, viscosity"
          extends Modelica.Icons.FunctionsPackage;
          constant Boolean excludeEnthalpyOfFormation = true "If true, enthalpy of formation Hf is not included in specific enthalpy h";
          constant Modelica.Media.Interfaces.Choices.ReferenceEnthalpy referenceChoice = Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt0K "Choice of reference enthalpy";
          constant Modelica.Media.Interfaces.Types.SpecificEnthalpy h_offset = 0.0 "User defined offset for reference enthalpy, if referenceChoice = UserDefined";
          constant Integer methodForThermalConductivity(min = 1, max = 2) = 1;

          function cp_T "Compute specific heat capacity at constant pressure from temperature and gas data"
            extends Modelica.Icons.Function;
            input IdealGases.Common.DataRecord data "Ideal gas data";
            input SI.Temperature T "Temperature";
            output SI.SpecificHeatCapacity cp "Specific heat capacity at temperature T";
          algorithm
            cp := smooth(0, if T < data.Tlimit then data.R * (1 / (T * T) * (data.alow[1] + T * (data.alow[2] + T * (1. * data.alow[3] + T * (data.alow[4] + T * (data.alow[5] + T * (data.alow[6] + data.alow[7] * T))))))) else data.R * (1 / (T * T) * (data.ahigh[1] + T * (data.ahigh[2] + T * (1. * data.ahigh[3] + T * (data.ahigh[4] + T * (data.ahigh[5] + T * (data.ahigh[6] + data.ahigh[7] * T))))))));
            annotation(Inline = true, smoothOrder = 2);
          end cp_T;

          function h_T "Compute specific enthalpy from temperature and gas data; reference is decided by the
              refChoice input, or by the referenceChoice package constant by default"
            import Modelica.Media.Interfaces.Choices;
            extends Modelica.Icons.Function;
            input IdealGases.Common.DataRecord data "Ideal gas data";
            input SI.Temperature T "Temperature";
            input Boolean exclEnthForm = excludeEnthalpyOfFormation "If true, enthalpy of formation Hf is not included in specific enthalpy h";
            input Modelica.Media.Interfaces.Choices.ReferenceEnthalpy refChoice = referenceChoice "Choice of reference enthalpy";
            input SI.SpecificEnthalpy h_off = h_offset "User defined offset for reference enthalpy, if referenceChoice = UserDefined";
            output SI.SpecificEnthalpy h "Specific enthalpy at temperature T";
          algorithm
            h := smooth(0, (if T < data.Tlimit then data.R * (((-data.alow[1]) + T * (data.blow[1] + data.alow[2] * Math.log(T) + T * (1. * data.alow[3] + T * (0.5 * data.alow[4] + T * (1 / 3 * data.alow[5] + T * (0.25 * data.alow[6] + 0.2 * data.alow[7] * T)))))) / T) else data.R * (((-data.ahigh[1]) + T * (data.bhigh[1] + data.ahigh[2] * Math.log(T) + T * (1. * data.ahigh[3] + T * (0.5 * data.ahigh[4] + T * (1 / 3 * data.ahigh[5] + T * (0.25 * data.ahigh[6] + 0.2 * data.ahigh[7] * T)))))) / T)) + (if exclEnthForm then -data.Hf else 0.0) + (if refChoice == Choices.ReferenceEnthalpy.ZeroAt0K then data.H0 else 0.0) + (if refChoice == Choices.ReferenceEnthalpy.UserDefined then h_off else 0.0));
            annotation(Inline = false, smoothOrder = 2);
          end h_T;

          function s0_T "Compute specific entropy from temperature and gas data"
            extends Modelica.Icons.Function;
            input IdealGases.Common.DataRecord data "Ideal gas data";
            input SI.Temperature T "Temperature";
            output SI.SpecificEntropy s "Specific entropy at temperature T";
          algorithm
            s := if T < data.Tlimit then data.R * (data.blow[2] - 0.5 * data.alow[1] / (T * T) - data.alow[2] / T + data.alow[3] * Math.log(T) + T * (data.alow[4] + T * (0.5 * data.alow[5] + T * (1 / 3 * data.alow[6] + 0.25 * data.alow[7] * T)))) else data.R * (data.bhigh[2] - 0.5 * data.ahigh[1] / (T * T) - data.ahigh[2] / T + data.ahigh[3] * Math.log(T) + T * (data.ahigh[4] + T * (0.5 * data.ahigh[5] + T * (1 / 3 * data.ahigh[6] + 0.25 * data.ahigh[7] * T))));
            annotation(Inline = true, smoothOrder = 2);
          end s0_T;

          function dynamicViscosityLowPressure "Dynamic viscosity of low pressure gases"
            extends Modelica.Icons.Function;
            input SI.Temp_K T "Gas temperature";
            input SI.Temp_K Tc "Critical temperature of gas";
            input SI.MolarMass M "Molar mass of gas";
            input SI.MolarVolume Vc "Critical molar volume of gas";
            input Real w "Acentric factor of gas";
            input Modelica.Media.Interfaces.Types.DipoleMoment mu "Dipole moment of gas molecule";
            input Real k = 0.0 "Special correction for highly polar substances";
            output SI.DynamicViscosity eta "Dynamic viscosity of gas";
          protected
            parameter Real Const1_SI = 40.785 * 10 ^ (-9.5) "Constant in formula for eta converted to SI units";
            parameter Real Const2_SI = 131.3 / 1000.0 "Constant in formula for mur converted to SI units";
            Real mur = Const2_SI * mu / sqrt(Vc * Tc) "Dimensionless dipole moment of gas molecule";
            Real Fc = 1 - 0.2756 * w + 0.059035 * mur ^ 4 + k "Factor to account for molecular shape and polarities of gas";
            Real Tstar "Dimensionless temperature defined by equation below";
            Real Ov "Viscosity collision integral for the gas";
          algorithm
            Tstar := 1.2593 * T / Tc;
            Ov := 1.16145 * Tstar ^ (-0.14874) + 0.52487 * Modelica.Math.exp(-0.7732 * Tstar) + 2.16178 * Modelica.Math.exp(-2.43787 * Tstar);
            eta := Const1_SI * Fc * sqrt(M * T) / (Vc ^ (2 / 3) * Ov);
            annotation(smoothOrder = 2);
          end dynamicViscosityLowPressure;

          function thermalConductivityEstimate "Thermal conductivity of polyatomic gases (Eucken and Modified Eucken correlation)"
            extends Modelica.Icons.Function;
            input Modelica.Media.Interfaces.Types.SpecificHeatCapacity Cp "Constant pressure heat capacity";
            input Modelica.Media.Interfaces.Types.DynamicViscosity eta "Dynamic viscosity";
            input Integer method(min = 1, max = 2) = 1 "1: Eucken Method, 2: Modified Eucken Method";
            input IdealGases.Common.DataRecord data "Ideal gas data";
            output Modelica.Media.Interfaces.Types.ThermalConductivity lambda "Thermal conductivity [W/(m.k)]";
          algorithm
            lambda := if method == 1 then eta * (Cp - data.R + 9 / 4 * data.R) else eta * (Cp - data.R) * (1.32 + 1.77 / (Cp / data.R - 1.0));
            annotation(smoothOrder = 2);
          end thermalConductivityEstimate;
        end Functions;

        partial package SingleGasNasa "Medium model of an ideal gas based on NASA source"
          extends Interfaces.PartialPureSubstance(ThermoStates = Modelica.Media.Interfaces.Choices.IndependentVariables.pT, redeclare record FluidConstants = Modelica.Media.Interfaces.Types.IdealGas.FluidConstants, mediumName = data.name, substanceNames = {data.name}, singleState = false, Temperature(min = 200, max = 6000, start = 500, nominal = 500), SpecificEnthalpy(start = if Functions.referenceChoice == ReferenceEnthalpy.ZeroAt0K then data.H0 else if Functions.referenceChoice == ReferenceEnthalpy.UserDefined then Functions.h_offset else 0, nominal = 1.0e5), Density(start = 10, nominal = 10), AbsolutePressure(start = 10e5, nominal = 10e5));

          redeclare record extends ThermodynamicState "Thermodynamic state variables for ideal gases"
            AbsolutePressure p "Absolute pressure of medium";
            Temperature T "Temperature of medium";
          end ThermodynamicState;

          import Modelica.Math;
          import Modelica.Media.Interfaces.Choices.ReferenceEnthalpy;
          constant IdealGases.Common.DataRecord data "Data record of ideal gas substance";
          constant FluidConstants[nS] fluidConstants "Constant data for the fluid";

          redeclare model extends BaseProperties(T(stateSelect = if preferredMediumStates then StateSelect.prefer else StateSelect.default), p(stateSelect = if preferredMediumStates then StateSelect.prefer else StateSelect.default)) "Base properties of ideal gas medium"
          equation
            assert(T >= 200 and T <= 6000, "
          Temperature T (= " + String(T) + " K) is not in the allowed range
          200 K <= T <= 6000 K required from medium model \"" + mediumName + "\".
            ");
            MM = data.MM;
            R = data.R;
            h = Modelica.Media.IdealGases.Common.Functions.h_T(data, T, Modelica.Media.IdealGases.Common.Functions.excludeEnthalpyOfFormation, Modelica.Media.IdealGases.Common.Functions.referenceChoice, Modelica.Media.IdealGases.Common.Functions.h_offset);
            u = h - R * T;
            d = p / (R * T);
            state.T = T;
            state.p = p;
          end BaseProperties;

          redeclare function setState_pTX "Return thermodynamic state as function of p, T and composition X"
            extends Modelica.Icons.Function;
            input AbsolutePressure p "Pressure";
            input Temperature T "Temperature";
            input MassFraction[:] X = reference_X "Mass fractions";
            output ThermodynamicState state;
          algorithm
            state := ThermodynamicState(p = p, T = T);
            annotation(Inline = true, smoothOrder = 2);
          end setState_pTX;

          redeclare function setState_phX "Return thermodynamic state as function of p, h and composition X"
            extends Modelica.Icons.Function;
            input AbsolutePressure p "Pressure";
            input SpecificEnthalpy h "Specific enthalpy";
            input MassFraction[:] X = reference_X "Mass fractions";
            output ThermodynamicState state;
          algorithm
            state := ThermodynamicState(p = p, T = T_h(h));
            annotation(Inline = true, smoothOrder = 2);
          end setState_phX;

          redeclare function setState_psX "Return thermodynamic state as function of p, s and composition X"
            extends Modelica.Icons.Function;
            input AbsolutePressure p "Pressure";
            input SpecificEntropy s "Specific entropy";
            input MassFraction[:] X = reference_X "Mass fractions";
            output ThermodynamicState state;
          algorithm
            state := ThermodynamicState(p = p, T = T_ps(p, s));
            annotation(Inline = true, smoothOrder = 2);
          end setState_psX;

          redeclare function setState_dTX "Return thermodynamic state as function of d, T and composition X"
            extends Modelica.Icons.Function;
            input Density d "Density";
            input Temperature T "Temperature";
            input MassFraction[:] X = reference_X "Mass fractions";
            output ThermodynamicState state;
          algorithm
            state := ThermodynamicState(p = d * data.R * T, T = T);
            annotation(Inline = true, smoothOrder = 2);
          end setState_dTX;

          redeclare function extends setSmoothState "Return thermodynamic state so that it smoothly approximates: if x > 0 then state_a else state_b"
          algorithm
            state := ThermodynamicState(p = Media.Common.smoothStep(x, state_a.p, state_b.p, x_small), T = Media.Common.smoothStep(x, state_a.T, state_b.T, x_small));
            annotation(Inline = true, smoothOrder = 2);
          end setSmoothState;

          redeclare function extends pressure "Return pressure of ideal gas"
          algorithm
            p := state.p;
            annotation(Inline = true, smoothOrder = 2);
          end pressure;

          redeclare function extends temperature "Return temperature of ideal gas"
          algorithm
            T := state.T;
            annotation(Inline = true, smoothOrder = 2);
          end temperature;

          redeclare function extends density "Return density of ideal gas"
          algorithm
            d := state.p / (data.R * state.T);
            annotation(Inline = true, smoothOrder = 2);
          end density;

          redeclare function extends specificEnthalpy "Return specific enthalpy"
            extends Modelica.Icons.Function;
          algorithm
            h := Modelica.Media.IdealGases.Common.Functions.h_T(data, state.T);
            annotation(Inline = true, smoothOrder = 2);
          end specificEnthalpy;

          redeclare function extends specificInternalEnergy "Return specific internal energy"
            extends Modelica.Icons.Function;
          algorithm
            u := Modelica.Media.IdealGases.Common.Functions.h_T(data, state.T) - data.R * state.T;
            annotation(Inline = true, smoothOrder = 2);
          end specificInternalEnergy;

          redeclare function extends specificEntropy "Return specific entropy"
            extends Modelica.Icons.Function;
          algorithm
            s := Modelica.Media.IdealGases.Common.Functions.s0_T(data, state.T) - data.R * Modelica.Math.log(state.p / reference_p);
            annotation(Inline = true, smoothOrder = 2);
          end specificEntropy;

          redeclare function extends specificGibbsEnergy "Return specific Gibbs energy"
            extends Modelica.Icons.Function;
          algorithm
            g := Modelica.Media.IdealGases.Common.Functions.h_T(data, state.T) - state.T * specificEntropy(state);
            annotation(Inline = true, smoothOrder = 2);
          end specificGibbsEnergy;

          redeclare function extends specificHelmholtzEnergy "Return specific Helmholtz energy"
            extends Modelica.Icons.Function;
          algorithm
            f := Modelica.Media.IdealGases.Common.Functions.h_T(data, state.T) - data.R * state.T - state.T * specificEntropy(state);
            annotation(Inline = true, smoothOrder = 2);
          end specificHelmholtzEnergy;

          redeclare function extends specificHeatCapacityCp "Return specific heat capacity at constant pressure"
          algorithm
            cp := Modelica.Media.IdealGases.Common.Functions.cp_T(data, state.T);
            annotation(Inline = true, smoothOrder = 2);
          end specificHeatCapacityCp;

          redeclare function extends specificHeatCapacityCv "Compute specific heat capacity at constant volume from temperature and gas data"
          algorithm
            cv := Modelica.Media.IdealGases.Common.Functions.cp_T(data, state.T) - data.R;
            annotation(Inline = true, smoothOrder = 2);
          end specificHeatCapacityCv;

          redeclare function extends isentropicExponent "Return isentropic exponent"
          algorithm
            gamma := specificHeatCapacityCp(state) / specificHeatCapacityCv(state);
            annotation(Inline = true, smoothOrder = 2);
          end isentropicExponent;

          redeclare function extends velocityOfSound "Return velocity of sound"
            extends Modelica.Icons.Function;
          algorithm
            a := sqrt(max(0, data.R * state.T * Modelica.Media.IdealGases.Common.Functions.cp_T(data, state.T) / specificHeatCapacityCv(state)));
            annotation(Inline = true, smoothOrder = 2);
          end velocityOfSound;

          function isentropicEnthalpyApproximation "Approximate method of calculating h_is from upstream properties and downstream pressure"
            extends Modelica.Icons.Function;
            input SI.Pressure p2 "Downstream pressure";
            input ThermodynamicState state "Properties at upstream location";
            input Boolean exclEnthForm = Functions.excludeEnthalpyOfFormation "If true, enthalpy of formation Hf is not included in specific enthalpy h";
            input ReferenceEnthalpy refChoice = Functions.referenceChoice "Choice of reference enthalpy";
            input SpecificEnthalpy h_off = Functions.h_offset "User defined offset for reference enthalpy, if referenceChoice = UserDefined";
            output SI.SpecificEnthalpy h_is "Isentropic enthalpy";
          protected
            IsentropicExponent gamma = isentropicExponent(state) "Isentropic exponent";
          algorithm
            h_is := Modelica.Media.IdealGases.Common.Functions.h_T(data, state.T, exclEnthForm, refChoice, h_off) + gamma / (gamma - 1.0) * state.p / density(state) * ((p2 / state.p) ^ ((gamma - 1) / gamma) - 1.0);
            annotation(Inline = true, smoothOrder = 2);
          end isentropicEnthalpyApproximation;

          redeclare function extends isentropicEnthalpy "Return isentropic enthalpy"
            input Boolean exclEnthForm = Functions.excludeEnthalpyOfFormation "If true, enthalpy of formation Hf is not included in specific enthalpy h";
            input ReferenceEnthalpy refChoice = Functions.referenceChoice "Choice of reference enthalpy";
            input SpecificEnthalpy h_off = Functions.h_offset "User defined offset for reference enthalpy, if referenceChoice = UserDefined";
          algorithm
            h_is := isentropicEnthalpyApproximation(p_downstream, refState, exclEnthForm, refChoice, h_off);
            annotation(Inline = true, smoothOrder = 2);
          end isentropicEnthalpy;

          redeclare function extends isobaricExpansionCoefficient "Returns overall the isobaric expansion coefficient beta"
          algorithm
            beta := 1 / state.T;
            annotation(Inline = true, smoothOrder = 2);
          end isobaricExpansionCoefficient;

          redeclare function extends isothermalCompressibility "Returns overall the isothermal compressibility factor"
          algorithm
            kappa := 1.0 / state.p;
            annotation(Inline = true, smoothOrder = 2);
          end isothermalCompressibility;

          redeclare function extends density_derp_T "Returns the partial derivative of density with respect to pressure at constant temperature"
          algorithm
            ddpT := 1 / (state.T * data.R);
            annotation(Inline = true, smoothOrder = 2);
          end density_derp_T;

          redeclare function extends density_derT_p "Returns the partial derivative of density with respect to temperature at constant pressure"
          algorithm
            ddTp := -state.p / (state.T * state.T * data.R);
            annotation(Inline = true, smoothOrder = 2);
          end density_derT_p;

          redeclare function extends density_derX "Returns the partial derivative of density with respect to mass fractions at constant pressure and temperature"
          algorithm
            dddX := fill(0, nX);
            annotation(Inline = true, smoothOrder = 2);
          end density_derX;

          redeclare replaceable function extends dynamicViscosity "Dynamic viscosity"
          algorithm
            assert(fluidConstants[1].hasCriticalData, "Failed to compute dynamicViscosity: For the species \"" + mediumName + "\" no critical data is available.");
            assert(fluidConstants[1].hasDipoleMoment, "Failed to compute dynamicViscosity: For the species \"" + mediumName + "\" no critical data is available.");
            eta := Modelica.Media.IdealGases.Common.Functions.dynamicViscosityLowPressure(state.T, fluidConstants[1].criticalTemperature, fluidConstants[1].molarMass, fluidConstants[1].criticalMolarVolume, fluidConstants[1].acentricFactor, fluidConstants[1].dipoleMoment);
            annotation(smoothOrder = 2);
          end dynamicViscosity;

          redeclare replaceable function extends thermalConductivity "Thermal conductivity of gas"
            input Integer method = Functions.methodForThermalConductivity "1: Eucken Method, 2: Modified Eucken Method";
          algorithm
            assert(fluidConstants[1].hasCriticalData, "Failed to compute thermalConductivity: For the species \"" + mediumName + "\" no critical data is available.");
            lambda := Modelica.Media.IdealGases.Common.Functions.thermalConductivityEstimate(specificHeatCapacityCp(state), dynamicViscosity(state), method = method, data = data);
            annotation(smoothOrder = 2);
          end thermalConductivity;

          redeclare function extends molarMass "Return the molar mass of the medium"
          algorithm
            MM := data.MM;
            annotation(Inline = true, smoothOrder = 2);
          end molarMass;

          function T_h "Compute temperature from specific enthalpy"
            extends Modelica.Icons.Function;
            input SpecificEnthalpy h "Specific enthalpy";
            output Temperature T "Temperature";

          protected
            package Internal "Solve h(data,T) for T with given h (use only indirectly via temperature_phX)"
              extends Modelica.Media.Common.OneNonLinearEquation;

              redeclare record extends f_nonlinear_Data "Data to be passed to non-linear function"
                extends Modelica.Media.IdealGases.Common.DataRecord;
              end f_nonlinear_Data;

              redeclare function extends f_nonlinear
              algorithm
                y := Modelica.Media.IdealGases.Common.Functions.h_T(f_nonlinear_data, x);
              end f_nonlinear;

              redeclare function extends solve end solve;
            end Internal;
          algorithm
            T := Internal.solve(h, 200, 6000, 1.0e5, {1}, data);
          end T_h;

          function T_ps "Compute temperature from pressure and specific entropy"
            extends Modelica.Icons.Function;
            input AbsolutePressure p "Pressure";
            input SpecificEntropy s "Specific entropy";
            output Temperature T "Temperature";

          protected
            package Internal "Solve h(data,T) for T with given h (use only indirectly via temperature_phX)"
              extends Modelica.Media.Common.OneNonLinearEquation;

              redeclare record extends f_nonlinear_Data "Data to be passed to non-linear function"
                extends Modelica.Media.IdealGases.Common.DataRecord;
              end f_nonlinear_Data;

              redeclare function extends f_nonlinear
              algorithm
                y := Modelica.Media.IdealGases.Common.Functions.s0_T(f_nonlinear_data, x) - data.R * Modelica.Math.log(p / reference_p);
              end f_nonlinear;

              redeclare function extends solve end solve;
            end Internal;
          algorithm
            T := Internal.solve(s, 200, 6000, p, {1}, data);
          end T_ps;
        end SingleGasNasa;

        package FluidData "Critical data, dipole moments and related data"
          extends Modelica.Icons.Package;
          import Modelica.Media.Interfaces.PartialMixtureMedium;
          import Modelica.Media.IdealGases.Common.SingleGasesData;
          constant Modelica.Media.Interfaces.Types.IdealGas.FluidConstants N2(chemicalFormula = "N2", iupacName = "unknown", structureFormula = "unknown", casRegistryNumber = "7727-37-9", meltingPoint = 63.15, normalBoilingPoint = 77.35, criticalTemperature = 126.20, criticalPressure = 33.98e5, criticalMolarVolume = 90.10e-6, acentricFactor = 0.037, dipoleMoment = 0.0, molarMass = SingleGasesData.N2.MM, hasDipoleMoment = true, hasIdealGasHeatCapacity = true, hasCriticalData = true, hasAcentricFactor = true);
        end FluidData;

        package SingleGasesData "Ideal gas data based on the NASA Glenn coefficients"
          extends Modelica.Icons.Package;
          constant IdealGases.Common.DataRecord Air(name = "Air", MM = 0.0289651159, Hf = -4333.833858403446, H0 = 298609.6803431054, Tlimit = 1000, alow = {10099.5016, -196.827561, 5.00915511, -0.00576101373, 1.06685993e-005, -7.94029797e-009, 2.18523191e-012}, blow = {-176.796731, -3.921504225}, ahigh = {241521.443, -1257.8746, 5.14455867, -0.000213854179, 7.06522784e-008, -1.07148349e-011, 6.57780015e-016}, bhigh = {6462.26319, -8.147411905}, R = 287.0512249529787);
          constant IdealGases.Common.DataRecord Ar(name = "Ar", MM = 0.039948, Hf = 0, H0 = 155137.3785921698, Tlimit = 1000, alow = {0, 0, 2.5, 0, 0, 0, 0}, blow = {-745.375, 4.37967491}, ahigh = {20.10538475, -0.05992661069999999, 2.500069401, -3.99214116e-008, 1.20527214e-011, -1.819015576e-015, 1.078576636e-019}, bhigh = {-744.993961, 4.37918011}, R = 208.1323720837088);
          constant IdealGases.Common.DataRecord CH4(name = "CH4", MM = 0.01604246, Hf = -4650159.63885838, H0 = 624355.7409524474, Tlimit = 1000, alow = {-176685.0998, 2786.18102, -12.0257785, 0.0391761929, -3.61905443e-005, 2.026853043e-008, -4.976705489999999e-012}, blow = {-23313.1436, 89.0432275}, ahigh = {3730042.76, -13835.01485, 20.49107091, -0.001961974759, 4.72731304e-007, -3.72881469e-011, 1.623737207e-015}, bhigh = {75320.6691, -121.9124889}, R = 518.2791167938085);
          constant IdealGases.Common.DataRecord CH3OH(name = "CH3OH", MM = 0.03204186, Hf = -6271171.523750494, H0 = 356885.5553329301, Tlimit = 1000, alow = {-241664.2886, 4032.14719, -20.46415436, 0.0690369807, -7.59893269e-005, 4.59820836e-008, -1.158706744e-011}, blow = {-44332.61169999999, 140.014219}, ahigh = {3411570.76, -13455.00201, 22.61407623, -0.002141029179, 3.73005054e-007, -3.49884639e-011, 1.366073444e-015}, bhigh = {56360.8156, -127.7814279}, R = 259.4878075117987);
          constant IdealGases.Common.DataRecord CO(name = "CO", MM = 0.0280101, Hf = -3946262.098314536, H0 = 309570.6191695138, Tlimit = 1000, alow = {14890.45326, -292.2285939, 5.72452717, -0.008176235030000001, 1.456903469e-005, -1.087746302e-008, 3.027941827e-012}, blow = {-13031.31878, -7.85924135}, ahigh = {461919.725, -1944.704863, 5.91671418, -0.0005664282830000001, 1.39881454e-007, -1.787680361e-011, 9.62093557e-016}, bhigh = {-2466.261084, -13.87413108}, R = 296.8383547363272);
          constant IdealGases.Common.DataRecord CO2(name = "CO2", MM = 0.0440095, Hf = -8941478.544405185, H0 = 212805.6215135368, Tlimit = 1000, alow = {49436.5054, -626.411601, 5.30172524, 0.002503813816, -2.127308728e-007, -7.68998878e-010, 2.849677801e-013}, blow = {-45281.9846, -7.04827944}, ahigh = {117696.2419, -1788.791477, 8.29152319, -9.22315678e-005, 4.86367688e-009, -1.891053312e-012, 6.330036589999999e-016}, bhigh = {-39083.5059, -26.52669281}, R = 188.9244822140674);
          constant IdealGases.Common.DataRecord C2H2_vinylidene(name = "C2H2_vinylidene", MM = 0.02603728, Hf = 15930556.80163212, H0 = 417638.4015534649, Tlimit = 1000, alow = {-14660.42239, 278.9475593, 1.276229776, 0.01395015463, -1.475702649e-005, 9.476298110000001e-009, -2.567602217e-012}, blow = {47361.1018, 16.58225704}, ahigh = {1940838.725, -6892.718150000001, 13.39582494, -0.0009368968669999999, 1.470804368e-007, -1.220040365e-011, 4.12239166e-016}, bhigh = {91071.1293, -63.3750293}, R = 319.3295152181795);
          constant IdealGases.Common.DataRecord C2H4(name = "C2H4", MM = 0.02805316, Hf = 1871446.924339362, H0 = 374955.5843263291, Tlimit = 1000, alow = {-116360.5836, 2554.85151, -16.09746428, 0.0662577932, -7.885081859999999e-005, 5.12522482e-008, -1.370340031e-011}, blow = {-6176.19107, 109.3338343}, ahigh = {3408763.67, -13748.47903, 23.65898074, -0.002423804419, 4.43139566e-007, -4.35268339e-011, 1.775410633e-015}, bhigh = {88204.2938, -137.1278108}, R = 296.3827247982046);
          constant IdealGases.Common.DataRecord C2H6(name = "C2H6", MM = 0.03006904, Hf = -2788633.890539904, H0 = 395476.3437741943, Tlimit = 1000, alow = {-186204.4161, 3406.19186, -19.51705092, 0.0756583559, -8.20417322e-005, 5.0611358e-008, -1.319281992e-011}, blow = {-27029.3289, 129.8140496}, ahigh = {5025782.13, -20330.22397, 33.2255293, -0.00383670341, 7.23840586e-007, -7.3191825e-011, 3.065468699e-015}, bhigh = {111596.395, -203.9410584}, R = 276.5127187299628);
          constant IdealGases.Common.DataRecord C2H5OH(name = "C2H5OH", MM = 0.04606844, Hf = -5100020.751733725, H0 = 315659.1801241805, Tlimit = 1000, alow = {-234279.1392, 4479.18055, -27.44817302, 0.1088679162, -0.0001305309334, 8.437346399999999e-008, -2.234559017e-011}, blow = {-50222.29, 176.4829211}, ahigh = {4694817.65, -19297.98213, 34.4758404, -0.00323616598, 5.78494772e-007, -5.56460027e-011, 2.2262264e-015}, bhigh = {86016.22709999999, -203.4801732}, R = 180.4808671619877);
          constant IdealGases.Common.DataRecord C3H6_propylene(name = "C3H6_propylene", MM = 0.04207974, Hf = 475288.1077687267, H0 = 322020.9535515191, Tlimit = 1000, alow = {-191246.2174, 3542.07424, -21.14878626, 0.0890148479, -0.0001001429154, 6.267959389999999e-008, -1.637870781e-011}, blow = {-15299.61824, 140.7641382}, ahigh = {5017620.34, -20860.84035, 36.4415634, -0.00388119117, 7.27867719e-007, -7.321204500000001e-011, 3.052176369e-015}, bhigh = {126124.5355, -219.5715757}, R = 197.588483198803);
          constant IdealGases.Common.DataRecord C3H8(name = "C3H8", MM = 0.04409562, Hf = -2373931.923397381, H0 = 334301.1845620949, Tlimit = 1000, alow = {-243314.4337, 4656.27081, -29.39466091, 0.1188952745, -0.0001376308269, 8.814823909999999e-008, -2.342987994e-011}, blow = {-35403.3527, 184.1749277}, ahigh = {6420731.680000001, -26597.91134, 45.3435684, -0.00502066392, 9.471216939999999e-007, -9.57540523e-011, 4.00967288e-015}, bhigh = {145558.2459, -281.8374734}, R = 188.5555073270316);
          constant IdealGases.Common.DataRecord C4H8_1_butene(name = "C4H8_1_butene", MM = 0.05610631999999999, Hf = -9624.584182316718, H0 = 305134.9651875226, Tlimit = 1000, alow = {-272149.2014, 5100.079250000001, -31.8378625, 0.1317754442, -0.0001527359339, 9.714761109999999e-008, -2.56020447e-011}, blow = {-25230.96386, 200.6932108}, ahigh = {6257948.609999999, -26603.76305, 47.6492005, -0.00438326711, 7.12883844e-007, -5.991020839999999e-011, 2.051753504e-015}, bhigh = {156925.2657, -291.3869761}, R = 148.1913623991023);
          constant IdealGases.Common.DataRecord C4H10_n_butane(name = "C4H10_n_butane", MM = 0.0581222, Hf = -2164233.28779709, H0 = 330832.0228759407, Tlimit = 1000, alow = {-317587.254, 6176.331819999999, -38.9156212, 0.1584654284, -0.0001860050159, 1.199676349e-007, -3.20167055e-011}, blow = {-45403.63390000001, 237.9488665}, ahigh = {7682322.45, -32560.5151, 57.3673275, -0.00619791681, 1.180186048e-006, -1.221893698e-010, 5.250635250000001e-015}, bhigh = {177452.656, -358.791876}, R = 143.0515706563069);
          constant IdealGases.Common.DataRecord C5H10_1_pentene(name = "C5H10_1_pentene", MM = 0.07013290000000001, Hf = -303423.9279995551, H0 = 309127.3852927798, Tlimit = 1000, alow = {-534054.813, 9298.917380000001, -56.6779245, 0.2123100266, -0.000257129829, 1.666834304e-007, -4.43408047e-011}, blow = {-47906.8218, 339.60364}, ahigh = {3744014.97, -21044.85321, 47.3612699, -0.00042442012, -3.89897505e-008, 1.367074243e-011, -9.31319423e-016}, bhigh = {115409.1373, -278.6177449000001}, R = 118.5530899192818);
          constant IdealGases.Common.DataRecord C5H12_n_pentane(name = "C5H12_n_pentane", MM = 0.07214878, Hf = -2034130.029641527, H0 = 335196.2430965569, Tlimit = 1000, alow = {-276889.4625, 5834.28347, -36.1754148, 0.1533339707, -0.0001528395882, 8.191092e-008, -1.792327902e-011}, blow = {-46653.7525, 226.5544053}, ahigh = {-2530779.286, -8972.59326, 45.3622326, -0.002626989916, 3.135136419e-006, -5.31872894e-010, 2.886896868e-014}, bhigh = {14846.16529, -251.6550384}, R = 115.2406457877736);
          constant IdealGases.Common.DataRecord C6H6(name = "C6H6", MM = 0.07811184, Hf = 1061042.730525872, H0 = 181735.4577743912, Tlimit = 1000, alow = {-167734.0902, 4404.50004, -37.1737791, 0.1640509559, -0.0002020812374, 1.307915264e-007, -3.4442841e-011}, blow = {-10354.55401, 216.9853345}, ahigh = {4538575.72, -22605.02547, 46.940073, -0.004206676830000001, 7.90799433e-007, -7.9683021e-011, 3.32821208e-015}, bhigh = {139146.4686, -286.8751333}, R = 106.4431717393932);
          constant IdealGases.Common.DataRecord C6H12_1_hexene(name = "C6H12_1_hexene", MM = 0.08415948000000001, Hf = -498458.4030224521, H0 = 311788.9986962847, Tlimit = 1000, alow = {-666883.165, 11768.64939, -72.70998330000001, 0.2709398396, -0.00033332464, 2.182347097e-007, -5.85946882e-011}, blow = {-62157.8054, 428.682564}, ahigh = {733290.696, -14488.48641, 46.7121549, 0.00317297847, -5.24264652e-007, 4.28035582e-011, -1.472353254e-015}, bhigh = {66977.4041, -262.3643854}, R = 98.79424159940152);
          constant IdealGases.Common.DataRecord C6H14_n_hexane(name = "C6H14_n_hexane", MM = 0.08617535999999999, Hf = -1936980.593988816, H0 = 333065.0431863586, Tlimit = 1000, alow = {-581592.67, 10790.97724, -66.3394703, 0.2523715155, -0.0002904344705, 1.802201514e-007, -4.617223680000001e-011}, blow = {-72715.4457, 393.828354}, ahigh = {-3106625.684, -7346.087920000001, 46.94131760000001, 0.001693963977, 2.068996667e-006, -4.21214168e-010, 2.452345845e-014}, bhigh = {523.750312, -254.9967718}, R = 96.48317105956971);
          constant IdealGases.Common.DataRecord C7H14_1_heptene(name = "C7H14_1_heptene", MM = 0.09818605999999999, Hf = -639194.6066478277, H0 = 313588.3036756949, Tlimit = 1000, alow = {-744940.284, 13321.79893, -82.81694379999999, 0.3108065994, -0.000378677992, 2.446841042e-007, -6.488763869999999e-011}, blow = {-72178.8501, 485.667149}, ahigh = {-1927608.174, -9125.024420000002, 47.4817797, 0.00606766053, -8.684859080000001e-007, 5.81399526e-011, -1.473979569e-015}, bhigh = {26009.14656, -256.2880707}, R = 84.68077851377274);
          constant IdealGases.Common.DataRecord C7H16_n_heptane(name = "C7H16_n_heptane", MM = 0.10020194, Hf = -1874015.612871368, H0 = 331540.487140269, Tlimit = 1000, alow = {-612743.289, 11840.85437, -74.87188599999999, 0.2918466052, -0.000341679549, 2.159285269e-007, -5.65585273e-011}, blow = {-80134.0894, 440.721332}, ahigh = {9135632.469999999, -39233.1969, 78.8978085, -0.00465425193, 2.071774142e-006, -3.4425393e-010, 1.976834775e-014}, bhigh = {205070.8295, -485.110402}, R = 82.97715593131233);
          constant IdealGases.Common.DataRecord C8H10_ethylbenz(name = "C8H10_ethylbenz", MM = 0.106165, Hf = 281825.4603682946, H0 = 209862.0072528611, Tlimit = 1000, alow = {-469494, 9307.16836, -65.2176947, 0.2612080237, -0.000318175348, 2.051355473e-007, -5.40181735e-011}, blow = {-40738.7021, 378.090436}, ahigh = {5551564.100000001, -28313.80598, 60.6124072, 0.001042112857, -1.327426719e-006, 2.166031743e-010, -1.142545514e-014}, bhigh = {164224.1062, -369.176982}, R = 78.31650732350586);
          constant IdealGases.Common.DataRecord C8H18_n_octane(name = "C8H18_n_octane", MM = 0.11422852, Hf = -1827477.060895125, H0 = 330740.51909278, Tlimit = 1000, alow = {-698664.715, 13385.01096, -84.1516592, 0.327193666, -0.000377720959, 2.339836988e-007, -6.01089265e-011}, blow = {-90262.2325, 493.922214}, ahigh = {6365406.949999999, -31053.64657, 69.6916234, 0.01048059637, -4.12962195e-006, 5.543226319999999e-010, -2.651436499e-014}, bhigh = {150096.8785, -416.989565}, R = 72.78805678301707);
          constant IdealGases.Common.DataRecord CL2(name = "CL2", MM = 0.07090600000000001, Hf = 0, H0 = 129482.8364313316, Tlimit = 1000, alow = {34628.1517, -554.7126520000001, 6.20758937, -0.002989632078, 3.17302729e-006, -1.793629562e-009, 4.260043590000001e-013}, blow = {1534.069331, -9.438331107}, ahigh = {6092569.42, -19496.27662, 28.54535795, -0.01449968764, 4.46389077e-006, -6.35852586e-010, 3.32736029e-014}, bhigh = {121211.7724, -169.0778824}, R = 117.2604857134798);
          constant IdealGases.Common.DataRecord F2(name = "F2", MM = 0.0379968064, Hf = 0, H0 = 232259.1511269747, Tlimit = 1000, alow = {10181.76308, 22.74241183, 1.97135304, 0.008151604010000001, -1.14896009e-005, 7.95865253e-009, -2.167079526e-012}, blow = {-958.6943, 11.30600296}, ahigh = {-2941167.79, 9456.5977, -7.73861615, 0.00764471299, -2.241007605e-006, 2.915845236e-010, -1.425033974e-014}, bhigh = {-60710.0561, 84.23835080000001}, R = 218.8202848542556);
          constant IdealGases.Common.DataRecord H2(name = "H2", MM = 0.00201588, Hf = 0, H0 = 4200697.462150524, Tlimit = 1000, alow = {40783.2321, -800.918604, 8.21470201, -0.01269714457, 1.753605076e-005, -1.20286027e-008, 3.36809349e-012}, blow = {2682.484665, -30.43788844}, ahigh = {560812.801, -837.150474, 2.975364532, 0.001252249124, -3.74071619e-007, 5.936625200000001e-011, -3.6069941e-015}, bhigh = {5339.82441, -2.202774769}, R = 4124.487568704486);
          constant IdealGases.Common.DataRecord H2O(name = "H2O", MM = 0.01801528, Hf = -13423382.81725291, H0 = 549760.6476280135, Tlimit = 1000, alow = {-39479.6083, 575.573102, 0.931782653, 0.00722271286, -7.34255737e-006, 4.95504349e-009, -1.336933246e-012}, blow = {-33039.7431, 17.24205775}, ahigh = {1034972.096, -2412.698562, 4.64611078, 0.002291998307, -6.836830479999999e-007, 9.426468930000001e-011, -4.82238053e-015}, bhigh = {-13842.86509, -7.97814851}, R = 461.5233290850878);
          constant IdealGases.Common.DataRecord He(name = "He", MM = 0.004002602, Hf = 0, H0 = 1548349.798456104, Tlimit = 1000, alow = {0, 0, 2.5, 0, 0, 0, 0}, blow = {-745.375, 0.9287239740000001}, ahigh = {0, 0, 2.5, 0, 0, 0, 0}, bhigh = {-745.375, 0.9287239740000001}, R = 2077.26673798694);
          constant IdealGases.Common.DataRecord NH3(name = "NH3", MM = 0.01703052, Hf = -2697510.117130892, H0 = 589713.1150428759, Tlimit = 1000, alow = {-76812.26149999999, 1270.951578, -3.89322913, 0.02145988418, -2.183766703e-005, 1.317385706e-008, -3.33232206e-012}, blow = {-12648.86413, 43.66014588}, ahigh = {2452389.535, -8040.89424, 12.71346201, -0.000398018658, 3.55250275e-008, 2.53092357e-012, -3.32270053e-016}, bhigh = {43861.91959999999, -64.62330602}, R = 488.2101075011215);
          constant IdealGases.Common.DataRecord NO(name = "NO", MM = 0.0300061, Hf = 3041758.509103149, H0 = 305908.1320131574, Tlimit = 1000, alow = {-11439.16503, 153.6467592, 3.43146873, -0.002668592368, 8.48139912e-006, -7.685111050000001e-009, 2.386797655e-012}, blow = {9098.214410000001, 6.72872549}, ahigh = {223901.8716, -1289.651623, 5.43393603, -0.00036560349, 9.880966450000001e-008, -1.416076856e-011, 9.380184619999999e-016}, bhigh = {17503.17656, -8.50166909}, R = 277.0927244793559);
          constant IdealGases.Common.DataRecord NO2(name = "NO2", MM = 0.0460055, Hf = 743237.6346306421, H0 = 221890.3174620426, Tlimit = 1000, alow = {-56420.3878, 963.308572, -2.434510974, 0.01927760886, -1.874559328e-005, 9.145497730000001e-009, -1.777647635e-012}, blow = {-1547.925037, 40.6785121}, ahigh = {721300.157, -3832.6152, 11.13963285, -0.002238062246, 6.54772343e-007, -7.6113359e-011, 3.32836105e-015}, bhigh = {25024.97403, -43.0513004}, R = 180.7277825477389);
          constant IdealGases.Common.DataRecord N2(name = "N2", MM = 0.0280134, Hf = 0, H0 = 309498.4543111511, Tlimit = 1000, alow = {22103.71497, -381.846182, 6.08273836, -0.00853091441, 1.384646189e-005, -9.62579362e-009, 2.519705809e-012}, blow = {710.846086, -10.76003744}, ahigh = {587712.406, -2239.249073, 6.06694922, -0.00061396855, 1.491806679e-007, -1.923105485e-011, 1.061954386e-015}, bhigh = {12832.10415, -15.86640027}, R = 296.8033869505308);
          constant IdealGases.Common.DataRecord N2O(name = "N2O", MM = 0.0440128, Hf = 1854006.107314236, H0 = 217685.1961247637, Tlimit = 1000, alow = {42882.2597, -644.011844, 6.03435143, 0.0002265394436, 3.47278285e-006, -3.62774864e-009, 1.137969552e-012}, blow = {11794.05506, -10.0312857}, ahigh = {343844.804, -2404.557558, 9.125636220000001, -0.000540166793, 1.315124031e-007, -1.4142151e-011, 6.38106687e-016}, bhigh = {21986.32638, -31.47805016}, R = 188.9103169986913);
          constant IdealGases.Common.DataRecord Ne(name = "Ne", MM = 0.0201797, Hf = 0, H0 = 307111.9986917546, Tlimit = 1000, alow = {0, 0, 2.5, 0, 0, 0, 0}, blow = {-745.375, 3.35532272}, ahigh = {0, 0, 2.5, 0, 0, 0, 0}, bhigh = {-745.375, 3.35532272}, R = 412.0215860493466);
          constant IdealGases.Common.DataRecord O2(name = "O2", MM = 0.0319988, Hf = 0, H0 = 271263.4223783392, Tlimit = 1000, alow = {-34255.6342, 484.700097, 1.119010961, 0.00429388924, -6.83630052e-007, -2.0233727e-009, 1.039040018e-012}, blow = {-3391.45487, 18.4969947}, ahigh = {-1037939.022, 2344.830282, 1.819732036, 0.001267847582, -2.188067988e-007, 2.053719572e-011, -8.193467050000001e-016}, bhigh = {-16890.10929, 17.38716506}, R = 259.8369938872708);
          constant IdealGases.Common.DataRecord SO2(name = "SO2", MM = 0.0640638, Hf = -4633037.690552231, H0 = 164650.3485587805, Tlimit = 1000, alow = {-53108.4214, 909.031167, -2.356891244, 0.02204449885, -2.510781471e-005, 1.446300484e-008, -3.36907094e-012}, blow = {-41137.52080000001, 40.45512519}, ahigh = {-112764.0116, -825.226138, 7.61617863, -0.000199932761, 5.65563143e-008, -5.45431661e-012, 2.918294102e-016}, bhigh = {-33513.0869, -16.55776085}, R = 129.7842463294403);
          constant IdealGases.Common.DataRecord SO3(name = "SO3", MM = 0.0800632, Hf = -4944843.573576874, H0 = 145990.9046852986, Tlimit = 1000, alow = {-39528.5529, 620.857257, -1.437731716, 0.02764126467, -3.144958662e-005, 1.792798e-008, -4.12638666e-012}, blow = {-51841.0617, 33.91331216}, ahigh = {-216692.3781, -1301.022399, 10.96287985, -0.000383710002, 8.466889039999999e-008, -9.70539929e-012, 4.49839754e-016}, bhigh = {-43982.83990000001, -36.55217314}, R = 103.8488594010732);
        end SingleGasesData;
      end Common;
    end IdealGases;

    package Incompressible "Medium model for T-dependent properties, defined by tables or polynomials"
      extends Modelica.Icons.VariantsPackage;
      import Modelica.Constants;
      import Modelica.Math;

      package Common "Common data structures"
        extends Modelica.Icons.Package;

        record BaseProps_Tpoly "Fluid state record"
          extends Modelica.Icons.Record;
          SI.Temperature T "Temperature";
          SI.Pressure p "Pressure";
        end BaseProps_Tpoly;
      end Common;

      package TableBased "Incompressible medium properties based on tables"
        import Poly = Modelica.Media.Incompressible.TableBased.Polynomials_Temp;
        extends Modelica.Media.Interfaces.PartialMedium(ThermoStates = if enthalpyOfT then Modelica.Media.Interfaces.Choices.IndependentVariables.T else Modelica.Media.Interfaces.Choices.IndependentVariables.pT, final reducedX = true, final fixedX = true, mediumName = "tableMedium", redeclare record ThermodynamicState = Common.BaseProps_Tpoly, singleState = true, reference_p = 1.013e5, Temperature(min = T_min, max = T_max));
        constant Boolean enthalpyOfT = true "True if enthalpy is approximated as a function of T only, (p-dependence neglected)";
        constant Boolean densityOfT = size(tableDensity, 1) > 1 "True if density is a function of temperature";
        constant Modelica.SIunits.Temperature T_min "Minimum temperature valid for medium model";
        constant Modelica.SIunits.Temperature T_max "Maximum temperature valid for medium model";
        constant Temperature T0 = 273.15 "Reference Temperature";
        constant SpecificEnthalpy h0 = 0 "Reference enthalpy at T0, reference_p";
        constant SpecificEntropy s0 = 0 "Reference entropy at T0, reference_p";
        constant MolarMass MM_const = 0.1 "Molar mass";
        constant Integer npol = 2 "Degree of polynomial used for fitting";
        constant Integer npolDensity = npol "Degree of polynomial used for fitting rho(T)";
        constant Integer npolHeatCapacity = npol "Degree of polynomial used for fitting Cp(T)";
        constant Integer npolViscosity = npol "Degree of polynomial used for fitting eta(T)";
        constant Integer npolVaporPressure = npol "Degree of polynomial used for fitting pVap(T)";
        constant Integer npolConductivity = npol "Degree of polynomial used for fitting lambda(T)";
        constant Integer neta = size(tableViscosity, 1) "Number of data points for viscosity";
        constant Real[:, 2] tableDensity "Table for rho(T)";
        constant Real[:, 2] tableHeatCapacity "Table for Cp(T)";
        constant Real[:, 2] tableViscosity "Table for eta(T)";
        constant Real[:, 2] tableVaporPressure "Table for pVap(T)";
        constant Real[:, 2] tableConductivity "Table for lambda(T)";
        constant Boolean TinK "True if T[K],Kelvin used for table temperatures";
        constant Boolean hasDensity = not size(tableDensity, 1) == 0 "True if table tableDensity is present";
        constant Boolean hasHeatCapacity = not size(tableHeatCapacity, 1) == 0 "True if table tableHeatCapacity is present";
        constant Boolean hasViscosity = not size(tableViscosity, 1) == 0 "True if table tableViscosity is present";
        constant Boolean hasVaporPressure = not size(tableVaporPressure, 1) == 0 "True if table tableVaporPressure is present";
        final constant Real[neta] invTK = if size(tableViscosity, 1) > 0 then if TinK then 1 ./ tableViscosity[:, 1] else 1 ./ Cv.from_degC(tableViscosity[:, 1]) else fill(0, neta);
        final constant Real[:] poly_rho = if hasDensity then Poly.fitting(tableDensity[:, 1], tableDensity[:, 2], npolDensity) else zeros(npolDensity + 1);
        final constant Real[:] poly_Cp = if hasHeatCapacity then Poly.fitting(tableHeatCapacity[:, 1], tableHeatCapacity[:, 2], npolHeatCapacity) else zeros(npolHeatCapacity + 1);
        final constant Real[:] poly_eta = if hasViscosity then Poly.fitting(invTK, Math.log(tableViscosity[:, 2]), npolViscosity) else zeros(npolViscosity + 1);
        final constant Real[:] poly_lam = if size(tableConductivity, 1) > 0 then Poly.fitting(tableConductivity[:, 1], tableConductivity[:, 2], npolConductivity) else zeros(npolConductivity + 1);

        redeclare model extends BaseProperties(final standardOrderComponents = true, p_bar = Cv.to_bar(p), T_degC(start = T_start - 273.15) = Cv.to_degC(T), T(start = T_start, stateSelect = if preferredMediumStates then StateSelect.prefer else StateSelect.default)) "Base properties of T dependent medium"
          SI.SpecificHeatCapacity cp "Specific heat capacity";
          parameter SI.Temperature T_start = 298.15 "Initial temperature";
        equation
          assert(hasDensity, "Medium " + mediumName + " can not be used without assigning tableDensity.");
          assert(T >= T_min and T <= T_max, "Temperature T (= " + String(T) + " K) is not in the allowed range (" + String(T_min) + " K <= T <= " + String(T_max) + " K) required from medium model \"" + mediumName + "\".");
          R = Modelica.Constants.R / MM_const;
          cp = Poly.evaluate(poly_Cp, if TinK then T else T_degC);
          h = specificEnthalpyOfT(p, T, densityOfT);
          u = h - (if singleState then reference_p / d else state.p / d);
          d = Poly.evaluate(poly_rho, if TinK then T else T_degC);
          state.T = T;
          state.p = p;
          MM = MM_const;
        end BaseProperties;

        redeclare function extends setState_pTX "Returns state record, given pressure and temperature"
        algorithm
          state := ThermodynamicState(p = p, T = T);
          annotation(smoothOrder = 3);
        end setState_pTX;

        redeclare function extends setState_dTX "Returns state record, given pressure and temperature"
        algorithm
          assert(false, "For incompressible media with d(T) only, state can not be set from density and temperature");
        end setState_dTX;

        redeclare function extends setState_phX "Returns state record, given pressure and specific enthalpy"
        algorithm
          state := ThermodynamicState(p = p, T = T_ph(p, h));
          annotation(Inline = true, smoothOrder = 3);
        end setState_phX;

        redeclare function extends setState_psX "Returns state record, given pressure and specific entropy"
        algorithm
          state := ThermodynamicState(p = p, T = T_ps(p, s));
          annotation(Inline = true, smoothOrder = 3);
        end setState_psX;

        redeclare function extends setSmoothState "Return thermodynamic state so that it smoothly approximates: if x > 0 then state_a else state_b"
        algorithm
          state := ThermodynamicState(p = Media.Common.smoothStep(x, state_a.p, state_b.p, x_small), T = Media.Common.smoothStep(x, state_a.T, state_b.T, x_small));
          annotation(Inline = true, smoothOrder = 3);
        end setSmoothState;

        redeclare function extends specificHeatCapacityCv "Specific heat capacity at constant volume (or pressure) of medium"
        algorithm
          assert(hasHeatCapacity, "Specific Heat Capacity, Cv, is not defined for medium " + mediumName + ".");
          cv := Poly.evaluate(poly_Cp, if TinK then state.T else state.T - 273.15);
          annotation(smoothOrder = 2);
        end specificHeatCapacityCv;

        redeclare function extends specificHeatCapacityCp "Specific heat capacity at constant volume (or pressure) of medium"
        algorithm
          assert(hasHeatCapacity, "Specific Heat Capacity, Cv, is not defined for medium " + mediumName + ".");
          cp := Poly.evaluate(poly_Cp, if TinK then state.T else state.T - 273.15);
          annotation(smoothOrder = 2);
        end specificHeatCapacityCp;

        redeclare function extends dynamicViscosity "Return dynamic viscosity as a function of the thermodynamic state record"
        algorithm
          assert(size(tableViscosity, 1) > 0, "DynamicViscosity, eta, is not defined for medium " + mediumName + ".");
          eta := Math.exp(Poly.evaluate(poly_eta, 1 / state.T));
          annotation(smoothOrder = 2);
        end dynamicViscosity;

        redeclare function extends thermalConductivity "Return thermal conductivity as a function of the thermodynamic state record"
        algorithm
          assert(size(tableConductivity, 1) > 0, "ThermalConductivity, lambda, is not defined for medium " + mediumName + ".");
          lambda := Poly.evaluate(poly_lam, if TinK then state.T else Cv.to_degC(state.T));
          annotation(smoothOrder = 2);
        end thermalConductivity;

        function s_T "Compute specific entropy"
          extends Modelica.Icons.Function;
          input Temperature T "Temperature";
          output SpecificEntropy s "Specific entropy";
        algorithm
          s := s0 + (if TinK then Poly.integralValue(poly_Cp[1:npol], T, T0) else Poly.integralValue(poly_Cp[1:npol], Cv.to_degC(T), Cv.to_degC(T0))) + Modelica.Math.log(T / T0) * Poly.evaluate(poly_Cp, if TinK then 0 else Modelica.Constants.T_zero);
          annotation(Inline = true, smoothOrder = 2);
        end s_T;

        redeclare function extends specificEntropy "Return specific entropy as a function of the thermodynamic state record"
        protected
          Integer npol = size(poly_Cp, 1) - 1;
        algorithm
          assert(hasHeatCapacity, "Specific Entropy, s(T), is not defined for medium " + mediumName + ".");
          s := s_T(state.T);
          annotation(smoothOrder = 2);
        end specificEntropy;

        function h_T "Compute specific enthalpy from temperature"
          import Modelica.SIunits.Conversions.to_degC;
          extends Modelica.Icons.Function;
          input SI.Temperature T "Temperature";
          output SI.SpecificEnthalpy h "Specific enthalpy at p, T";
        algorithm
          h := h0 + Poly.integralValue(poly_Cp, if TinK then T else Cv.to_degC(T), if TinK then T0 else Cv.to_degC(T0));
          annotation(derivative = h_T_der);
        end h_T;

        function h_T_der "Compute specific enthalpy from temperature"
          import Modelica.SIunits.Conversions.to_degC;
          extends Modelica.Icons.Function;
          input SI.Temperature T "Temperature";
          input Real dT "Temperature derivative";
          output Real dh "Derivative of Specific enthalpy at T";
        algorithm
          dh := Poly.evaluate(poly_Cp, if TinK then T else Cv.to_degC(T)) * dT;
          annotation(smoothOrder = 1);
        end h_T_der;

        function h_pT "Compute specific enthalpy from pressure and temperature"
          import Modelica.SIunits.Conversions.to_degC;
          extends Modelica.Icons.Function;
          input SI.Pressure p "Pressure";
          input SI.Temperature T "Temperature";
          input Boolean densityOfT = false "Include or neglect density derivative dependence of enthalpy";
          output SI.SpecificEnthalpy h "Specific enthalpy at p, T";
        algorithm
          h := h0 + Poly.integralValue(poly_Cp, if TinK then T else Cv.to_degC(T), if TinK then T0 else Cv.to_degC(T0)) + (p - reference_p) / Poly.evaluate(poly_rho, if TinK then T else Cv.to_degC(T)) * (if densityOfT then 1 + T / Poly.evaluate(poly_rho, if TinK then T else Cv.to_degC(T)) * Poly.derivativeValue(poly_rho, if TinK then T else Cv.to_degC(T)) else 1.0);
          annotation(smoothOrder = 2);
        end h_pT;

        redeclare function extends temperature "Return temperature as a function of the thermodynamic state record"
        algorithm
          T := state.T;
          annotation(Inline = true, smoothOrder = 2);
        end temperature;

        redeclare function extends pressure "Return pressure as a function of the thermodynamic state record"
        algorithm
          p := state.p;
          annotation(Inline = true, smoothOrder = 2);
        end pressure;

        redeclare function extends density "Return density as a function of the thermodynamic state record"
        algorithm
          d := Poly.evaluate(poly_rho, if TinK then state.T else Cv.to_degC(state.T));
          annotation(Inline = true, smoothOrder = 2);
        end density;

        redeclare function extends specificEnthalpy "Return specific enthalpy as a function of the thermodynamic state record"
        algorithm
          h := specificEnthalpyOfT(state.p, state.T);
          annotation(Inline = true, smoothOrder = 2);
        end specificEnthalpy;

        redeclare function extends specificInternalEnergy "Return specific internal energy as a function of the thermodynamic state record"
        algorithm
          u := specificEnthalpyOfT(state.p, state.T) - (if singleState then reference_p else state.p) / density(state);
          annotation(Inline = true, smoothOrder = 2);
        end specificInternalEnergy;

        function T_ph "Compute temperature from pressure and specific enthalpy"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SpecificEnthalpy h "Specific enthalpy";
          output Temperature T "Temperature";

        protected
          package Internal "Solve h(T) for T with given h (use only indirectly via temperature_phX)"
            extends Modelica.Media.Common.OneNonLinearEquation;

            redeclare record extends f_nonlinear_Data "Superfluous record, fix later when better structure of inverse functions exists"
              constant Real[5] dummy = {1, 2, 3, 4, 5};
            end f_nonlinear_Data;

            redeclare function extends f_nonlinear "P is smuggled in via vector"
            algorithm
              y := specificEnthalpyOfT(p, x);
            end f_nonlinear;
          end Internal;
        algorithm
          T := Internal.solve(h, T_min, T_max, p, {1}, Internal.f_nonlinear_Data());
          annotation(Inline = false, LateInline = true, inverse(h = specificEnthalpyOfT(p, T)));
        end T_ph;

        function T_ps "Compute temperature from pressure and specific enthalpy"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SpecificEntropy s "Specific entropy";
          output Temperature T "Temperature";

        protected
          package Internal "Solve h(T) for T with given h (use only indirectly via temperature_phX)"
            extends Modelica.Media.Common.OneNonLinearEquation;

            redeclare record extends f_nonlinear_Data "Superfluous record, fix later when better structure of inverse functions exists"
              constant Real[5] dummy = {1, 2, 3, 4, 5};
            end f_nonlinear_Data;

            redeclare function extends f_nonlinear "P is smuggled in via vector"
            algorithm
              y := s_T(x);
            end f_nonlinear;
          end Internal;
        algorithm
          T := Internal.solve(s, T_min, T_max, p, {1}, Internal.f_nonlinear_Data());
        end T_ps;

        package Polynomials_Temp "Temporary Functions operating on polynomials (including polynomial fitting); only to be used in Modelica.Media.Incompressible.TableBased"
          extends Modelica.Icons.Package;

          function evaluate "Evaluate polynomial at a given abscissa value"
            extends Modelica.Icons.Function;
            input Real[:] p "Polynomial coefficients (p[1] is coefficient of highest power)";
            input Real u "Abscissa value";
            output Real y "Value of polynomial at u";
          algorithm
            y := p[1];
            for j in 2:size(p, 1) loop
              y := p[j] + u * y;
            end for;
            annotation(derivative(zeroDerivative = p) = evaluate_der);
          end evaluate;

          function evaluateWithRange "Evaluate polynomial at a given abscissa value with linear extrapolation outside of the defined range"
            extends Modelica.Icons.Function;
            input Real[:] p "Polynomial coefficients (p[1] is coefficient of highest power)";
            input Real uMin "Polynomial valid in the range uMin .. uMax";
            input Real uMax "Polynomial valid in the range uMin .. uMax";
            input Real u "Abscissa value";
            output Real y "Value of polynomial at u. Outside of uMin,uMax, linear extrapolation is used";
          algorithm
            if u < uMin then
              y := evaluate(p, uMin) - evaluate_der(p, uMin, uMin - u);
            elseif u > uMax then
              y := evaluate(p, uMax) + evaluate_der(p, uMax, u - uMax);
            else
              y := evaluate(p, u);
            end if;
            annotation(derivative(zeroDerivative = p, zeroDerivative = uMin, zeroDerivative = uMax) = evaluateWithRange_der);
          end evaluateWithRange;

          function derivativeValue "Value of derivative of polynomial at abscissa value u"
            extends Modelica.Icons.Function;
            input Real[:] p "Polynomial coefficients (p[1] is coefficient of highest power)";
            input Real u "Abscissa value";
            output Real y "Value of derivative of polynomial at u";
          protected
            Integer n = size(p, 1);
          algorithm
            y := p[1] * (n - 1);
            for j in 2:size(p, 1) - 1 loop
              y := p[j] * (n - j) + u * y;
            end for;
            annotation(derivative(zeroDerivative = p) = derivativeValue_der);
          end derivativeValue;

          function secondDerivativeValue "Value of 2nd derivative of polynomial at abscissa value u"
            extends Modelica.Icons.Function;
            input Real[:] p "Polynomial coefficients (p[1] is coefficient of highest power)";
            input Real u "Abscissa value";
            output Real y "Value of 2nd derivative of polynomial at u";
          protected
            Integer n = size(p, 1);
          algorithm
            y := p[1] * (n - 1) * (n - 2);
            for j in 2:size(p, 1) - 2 loop
              y := p[j] * (n - j) * (n - j - 1) + u * y;
            end for;
          end secondDerivativeValue;

          function integralValue "Integral of polynomial p(u) from u_low to u_high"
            extends Modelica.Icons.Function;
            input Real[:] p "Polynomial coefficients";
            input Real u_high "High integrand value";
            input Real u_low = 0 "Low integrand value, default 0";
            output Real integral = 0.0 "Integral of polynomial p from u_low to u_high";
          protected
            Integer n = size(p, 1) "Degree of integrated polynomial";
            Real y_low = 0 "Value at lower integrand";
          algorithm
            for j in 1:n loop
              integral := u_high * (p[j] / (n - j + 1) + integral);
              y_low := u_low * (p[j] / (n - j + 1) + y_low);
            end for;
            integral := integral - y_low;
            annotation(derivative(zeroDerivative = p) = integralValue_der);
          end integralValue;

          function fitting "Computes the coefficients of a polynomial that fits a set of data points in a least-squares sense"
            extends Modelica.Icons.Function;
            input Real[:] u "Abscissa data values";
            input Real[size(u, 1)] y "Ordinate data values";
            input Integer n(min = 1) "Order of desired polynomial that fits the data points (u,y)";
            output Real[n + 1] p "Polynomial coefficients of polynomial that fits the date points";
          protected
            Real[size(u, 1), n + 1] V "Vandermonde matrix";
          algorithm
            V[:, n + 1] := ones(size(u, 1));
            for j in n:(-1):1 loop
              V[:, j] := {u[i] * V[i, j + 1] for i in 1:size(u, 1)};
            end for;
            p := Modelica.Math.Matrices.leastSquares(V, y);
          end fitting;

          function evaluate_der "Evaluate derivative of polynomial at a given abscissa value"
            extends Modelica.Icons.Function;
            input Real[:] p "Polynomial coefficients (p[1] is coefficient of highest power)";
            input Real u "Abscissa value";
            input Real du "Delta of abscissa value";
            output Real dy "Value of derivative of polynomial at u";
          protected
            Integer n = size(p, 1);
          algorithm
            dy := p[1] * (n - 1);
            for j in 2:size(p, 1) - 1 loop
              dy := p[j] * (n - j) + u * dy;
            end for;
            dy := dy * du;
          end evaluate_der;

          function evaluateWithRange_der "Evaluate derivative of polynomial at a given abscissa value with extrapolation outside of the defined range"
            extends Modelica.Icons.Function;
            input Real[:] p "Polynomial coefficients (p[1] is coefficient of highest power)";
            input Real uMin "Polynomial valid in the range uMin .. uMax";
            input Real uMax "Polynomial valid in the range uMin .. uMax";
            input Real u "Abscissa value";
            input Real du "Delta of abscissa value";
            output Real dy "Value of derivative of polynomial at u";
          algorithm
            if u < uMin then
              dy := evaluate_der(p, uMin, du);
            elseif u > uMax then
              dy := evaluate_der(p, uMax, du);
            else
              dy := evaluate_der(p, u, du);
            end if;
          end evaluateWithRange_der;

          function integralValue_der "Time derivative of integral of polynomial p(u) from u_low to u_high, assuming only u_high as time-dependent (Leibniz rule)"
            extends Modelica.Icons.Function;
            input Real[:] p "Polynomial coefficients";
            input Real u_high "High integrand value";
            input Real u_low = 0 "Low integrand value, default 0";
            input Real du_high "High integrand value";
            input Real du_low = 0 "Low integrand value, default 0";
            output Real dintegral = 0.0 "Integral of polynomial p from u_low to u_high";
          algorithm
            dintegral := evaluate(p, u_high) * du_high;
          end integralValue_der;

          function derivativeValue_der "Time derivative of derivative of polynomial"
            extends Modelica.Icons.Function;
            input Real[:] p "Polynomial coefficients (p[1] is coefficient of highest power)";
            input Real u "Abscissa value";
            input Real du "Delta of abscissa value";
            output Real dy "Time-derivative of derivative of polynomial w.r.t. input variable at u";
          protected
            Integer n = size(p, 1);
          algorithm
            dy := secondDerivativeValue(p, u) * du;
          end derivativeValue_der;
        end Polynomials_Temp;

      protected
        function specificEnthalpyOfT "Return specific enthalpy from pressure and temperature, taking the flag enthalpyOfT into account"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input Temperature T "Temperature";
          input Boolean densityOfT = false "Include or neglect density derivative dependence of enthalpy";
          output SpecificEnthalpy h "Specific enthalpy";
        algorithm
          h := if enthalpyOfT then h_T(T) else h_pT(p, T, densityOfT);
          annotation(Inline = true, smoothOrder = 2);
        end specificEnthalpyOfT;
      end TableBased;
    end Incompressible;
  end Media;

  package Math "Library of mathematical functions (e.g., sin, cos) and of functions operating on vectors and matrices"
    import SI = Modelica.SIunits;
    extends Modelica.Icons.Package;

    package Matrices "Library of functions operating on matrices"
      extends Modelica.Icons.Package;

      function leastSquares "Solve linear equation A*x = b (exactly if possible, or otherwise in a least square sense; A may be non-square and may be rank deficient)"
        extends Modelica.Icons.Function;
        input Real[:, :] A "Matrix A";
        input Real[size(A, 1)] b "Vector b";
        input Real rcond = 100 * Modelica.Constants.eps "Reciprocal condition number to estimate the rank of A";
        output Real[size(A, 2)] x "Vector x such that min|A*x-b|^2 if size(A,1) >= size(A,2) or min|x|^2 and A*x=b, if size(A,1) < size(A,2)";
        output Integer rank "Rank of A";
      protected
        Integer info;
        Real[max(size(A, 1), size(A, 2))] xx;
      algorithm
        if min(size(A)) > 0 then
          (xx, info, rank) := LAPACK.dgelsy_vec(A, b, rcond);
          x := xx[1:size(A, 2)];
          assert(info == 0, "Solving an overdetermined or underdetermined linear system\n" + "of equations with function \"Matrices.leastSquares\" failed.");
        else
          x := fill(0.0, size(A, 2));
        end if;
      end leastSquares;

      package LAPACK "Interface to LAPACK library (should usually not directly be used but only indirectly via Modelica.Math.Matrices)"
        extends Modelica.Icons.Package;

        function dgelsy_vec "Computes the minimum-norm solution to a real linear least squares problem with rank deficient A"
          extends Modelica.Icons.Function;
          input Real[:, :] A;
          input Real[size(A, 1)] b;
          input Real rcond = 0.0 "Reciprocal condition number to estimate rank";
          output Real[max(size(A, 1), size(A, 2))] x = cat(1, b, zeros(max(nrow, ncol) - nrow)) "solution is in first size(A,2) rows";
          output Integer info;
          output Integer rank "Effective rank of A";
        protected
          Integer nrow = size(A, 1);
          Integer ncol = size(A, 2);
          Integer nrhs = 1;
          Integer nx = max(nrow, ncol);
          Integer lwork = max(min(nrow, ncol) + 3 * ncol + 1, 2 * min(nrow, ncol) + 1);
          Real[max(min(size(A, 1), size(A, 2)) + 3 * size(A, 2) + 1, 2 * min(size(A, 1), size(A, 2)) + 1)] work;
          Real[size(A, 1), size(A, 2)] Awork = A;
          Integer[size(A, 2)] jpvt = zeros(ncol);
          external "FORTRAN 77" dgelsy(nrow, ncol, nrhs, Awork, nrow, x, nx, jpvt, rcond, rank, work, lwork, info) annotation(Library = "lapack");
        end dgelsy_vec;
      end LAPACK;
    end Matrices;

    package Icons "Icons for Math"
      extends Modelica.Icons.IconsPackage;

      partial function AxisLeft "Basic icon for mathematical function with y-axis on left side" end AxisLeft;

      partial function AxisCenter "Basic icon for mathematical function with y-axis in the center" end AxisCenter;
    end Icons;

    function asin "Inverse sine (-1 <= u <= 1)"
      extends Modelica.Math.Icons.AxisCenter;
      input Real u;
      output SI.Angle y;
      external "builtin" y = asin(u);
    end asin;

    function exp "Exponential, base e"
      extends Modelica.Math.Icons.AxisCenter;
      input Real u;
      output Real y;
      external "builtin" y = exp(u);
    end exp;

    function log "Natural (base e) logarithm (u shall be > 0)"
      extends Modelica.Math.Icons.AxisLeft;
      input Real u;
      output Real y;
      external "builtin" y = log(u);
    end log;
  end Math;

  package Utilities "Library of utility functions dedicated to scripting (operating on files, streams, strings, system)"
    extends Modelica.Icons.UtilitiesPackage;

    package Streams "Read from files and write to files"
      extends Modelica.Icons.FunctionsPackage;

      function error "Print error message and cancel all actions - in case of an unrecoverable error"
        extends Modelica.Icons.Function;
        input String string "String to be printed to error message window";
        external "C" ModelicaError(string) annotation(Library = "ModelicaExternalC");
      end error;
    end Streams;

    package Strings "Operations on strings"
      extends Modelica.Icons.FunctionsPackage;

      function length "Return length of string"
        extends Modelica.Icons.Function;
        input String string;
        output Integer result "Number of characters of string";
        external "C" result = ModelicaStrings_length(string) annotation(Library = "ModelicaExternalC");
      end length;

      function isEmpty "Return true if a string is empty (has only white space characters)"
        extends Modelica.Icons.Function;
        input String string;
        output Boolean result "True, if string is empty";
      protected
        Integer nextIndex;
        Integer len;
      algorithm
        nextIndex := Strings.Advanced.skipWhiteSpace(string);
        len := Strings.length(string);
        if len < 1 or nextIndex > len then
          result := true;
        else
          result := false;
        end if;
      end isEmpty;

      package Advanced "Advanced scanning functions"
        extends Modelica.Icons.FunctionsPackage;

        function skipWhiteSpace "Scan white space"
          extends Modelica.Icons.Function;
          input String string;
          input Integer startIndex(min = 1) = 1;
          output Integer nextIndex;
          external "C" nextIndex = ModelicaStrings_skipWhiteSpace(string, startIndex) annotation(Library = "ModelicaExternalC");
        end skipWhiteSpace;
      end Advanced;
    end Strings;
  end Utilities;

  package Constants "Library of mathematical constants and constants of nature (e.g., pi, eps, R, sigma)"
    import SI = Modelica.SIunits;
    import NonSI = Modelica.SIunits.Conversions.NonSIunits;
    extends Modelica.Icons.Package;
    final constant Real pi = 2 * Modelica.Math.asin(1.0);
    final constant Real eps = ModelicaServices.Machine.eps "Biggest number such that 1.0 + eps = 1.0";
    final constant Real inf = ModelicaServices.Machine.inf "Biggest Real number such that inf and -inf are representable on the machine";
    final constant SI.Velocity c = 299792458 "Speed of light in vacuum";
    final constant SI.FaradayConstant F = 9.648533289e4 "Faraday constant, C/mol (previous value: 9.64853399e4)";
    final constant Real R(final unit = "J/(mol.K)") = 8.3144598 "Molar gas constant (previous value: 8.314472)";
    final constant Real N_A(final unit = "1/mol") = 6.022140857e23 "Avogadro constant (previous value: 6.0221415e23)";
    final constant Real mue_0(final unit = "N/A2") = 4 * pi * 1.e-7 "Magnetic constant";
    final constant NonSI.Temperature_degC T_zero = -273.15 "Absolute zero temperature";
  end Constants;

  package Icons "Library of icons"
    partial class Information "Icon for general information packages" end Information;

    extends Icons.Package;

    partial package ExamplesPackage "Icon for packages containing runnable examples"
      extends Modelica.Icons.Package;
    end ExamplesPackage;

    partial package Package "Icon for standard packages" end Package;

    partial package BasesPackage "Icon for packages containing base classes"
      extends Modelica.Icons.Package;
    end BasesPackage;

    partial package VariantsPackage "Icon for package containing variants"
      extends Modelica.Icons.Package;
    end VariantsPackage;

    partial package InterfacesPackage "Icon for packages containing interfaces"
      extends Modelica.Icons.Package;
    end InterfacesPackage;

    partial package UtilitiesPackage "Icon for utility packages"
      extends Modelica.Icons.Package;
    end UtilitiesPackage;

    partial package TypesPackage "Icon for packages containing type definitions"
      extends Modelica.Icons.Package;
    end TypesPackage;

    partial package FunctionsPackage "Icon for packages containing functions"
      extends Modelica.Icons.Package;
    end FunctionsPackage;

    partial package IconsPackage "Icon for packages containing icons"
      extends Modelica.Icons.Package;
    end IconsPackage;

    partial package InternalPackage "Icon for an internal package (indicating that the package should not be directly utilized by user)" end InternalPackage;

    partial package MaterialPropertiesPackage "Icon for package containing property classes"
      extends Modelica.Icons.Package;
    end MaterialPropertiesPackage;

    partial class MaterialProperty "Icon for property classes" end MaterialProperty;

    partial function Function "Icon for functions" end Function;

    partial record Record "Icon for records" end Record;
  end Icons;

  package SIunits "Library of type and unit definitions based on SI units according to ISO 31-1992"
    extends Modelica.Icons.Package;

    package Icons "Icons for SIunits"
      extends Modelica.Icons.IconsPackage;

      partial function Conversion "Base icon for conversion functions" end Conversion;
    end Icons;

    package Conversions "Conversion functions to/from non SI units and type definitions of non SI units"
      extends Modelica.Icons.Package;

      package NonSIunits "Type definitions of non SI units"
        extends Modelica.Icons.Package;
        type Temperature_degC = Real(final quantity = "ThermodynamicTemperature", final unit = "degC") "Absolute temperature in degree Celsius (for relative temperature use SIunits.TemperatureDifference)" annotation(absoluteValue = true);
        type Pressure_bar = Real(final quantity = "Pressure", final unit = "bar") "Absolute pressure in bar";
      end NonSIunits;

      function to_degC "Convert from Kelvin to degCelsius"
        extends Modelica.SIunits.Icons.Conversion;
        input Temperature Kelvin "Kelvin value";
        output NonSIunits.Temperature_degC Celsius "Celsius value";
      algorithm
        Celsius := Kelvin + Modelica.Constants.T_zero;
        annotation(Inline = true);
      end to_degC;

      function from_degC "Convert from degCelsius to Kelvin"
        extends Modelica.SIunits.Icons.Conversion;
        input NonSIunits.Temperature_degC Celsius "Celsius value";
        output Temperature Kelvin "Kelvin value";
      algorithm
        Kelvin := Celsius - Modelica.Constants.T_zero;
        annotation(Inline = true);
      end from_degC;

      function to_bar "Convert from Pascal to bar"
        extends Modelica.SIunits.Icons.Conversion;
        input Pressure Pa "Pascal value";
        output NonSIunits.Pressure_bar bar "bar value";
      algorithm
        bar := Pa / 1e5;
        annotation(Inline = true);
      end to_bar;
    end Conversions;

    type Angle = Real(final quantity = "Angle", final unit = "rad", displayUnit = "deg");
    type Length = Real(final quantity = "Length", final unit = "m");
    type Position = Length;
    type Distance = Length(min = 0);
    type Radius = Length(min = 0);
    type Diameter = Length(min = 0);
    type Area = Real(final quantity = "Area", final unit = "m2");
    type Volume = Real(final quantity = "Volume", final unit = "m3");
    type Time = Real(final quantity = "Time", final unit = "s");
    type AngularVelocity = Real(final quantity = "AngularVelocity", final unit = "rad/s");
    type Velocity = Real(final quantity = "Velocity", final unit = "m/s");
    type Acceleration = Real(final quantity = "Acceleration", final unit = "m/s2");
    type Frequency = Real(final quantity = "Frequency", final unit = "Hz");
    type Mass = Real(quantity = "Mass", final unit = "kg", min = 0);
    type Density = Real(final quantity = "Density", final unit = "kg/m3", displayUnit = "g/cm3", min = 0.0);
    type SpecificVolume = Real(final quantity = "SpecificVolume", final unit = "m3/kg", min = 0.0);
    type MomentOfInertia = Real(final quantity = "MomentOfInertia", final unit = "kg.m2");
    type Inertia = MomentOfInertia;
    type Torque = Real(final quantity = "Torque", final unit = "N.m");
    type Pressure = Real(final quantity = "Pressure", final unit = "Pa", displayUnit = "bar");
    type AbsolutePressure = Pressure(min = 0.0, nominal = 1e5);
    type PressureDifference = Pressure;
    type DynamicViscosity = Real(final quantity = "DynamicViscosity", final unit = "Pa.s", min = 0);
    type Energy = Real(final quantity = "Energy", final unit = "J");
    type Power = Real(final quantity = "Power", final unit = "W");
    type MassFlowRate = Real(quantity = "MassFlowRate", final unit = "kg/s");
    type MomentumFlux = Real(final quantity = "MomentumFlux", final unit = "N");
    type ThermodynamicTemperature = Real(final quantity = "ThermodynamicTemperature", final unit = "K", min = 0.0, start = 288.15, nominal = 300, displayUnit = "degC") "Absolute temperature (use type TemperatureDifference for relative temperatures)" annotation(absoluteValue = true);
    type Temp_K = ThermodynamicTemperature;
    type Temperature = ThermodynamicTemperature;
    type Compressibility = Real(final quantity = "Compressibility", final unit = "1/Pa");
    type IsothermalCompressibility = Compressibility;
    type ThermalConductivity = Real(final quantity = "ThermalConductivity", final unit = "W/(m.K)");
    type HeatCapacity = Real(final quantity = "HeatCapacity", final unit = "J/K");
    type SpecificHeatCapacity = Real(final quantity = "SpecificHeatCapacity", final unit = "J/(kg.K)");
    type RatioOfSpecificHeatCapacities = Real(final quantity = "RatioOfSpecificHeatCapacities", final unit = "1");
    type Entropy = Real(final quantity = "Entropy", final unit = "J/K");
    type SpecificEntropy = Real(final quantity = "SpecificEntropy", final unit = "J/(kg.K)");
    type SpecificEnergy = Real(final quantity = "SpecificEnergy", final unit = "J/kg");
    type SpecificEnthalpy = SpecificEnergy;
    type DerDensityByPressure = Real(final unit = "s2/m2");
    type DerDensityByTemperature = Real(final unit = "kg/(m3.K)");
    type ElectricCharge = Real(final quantity = "ElectricCharge", final unit = "C");
    type AmountOfSubstance = Real(final quantity = "AmountOfSubstance", final unit = "mol", min = 0);
    type MolarMass = Real(final quantity = "MolarMass", final unit = "kg/mol", min = 0);
    type MolarVolume = Real(final quantity = "MolarVolume", final unit = "m3/mol", min = 0);
    type MassFraction = Real(final quantity = "MassFraction", final unit = "1", min = 0, max = 1);
    type MoleFraction = Real(final quantity = "MoleFraction", final unit = "1", min = 0, max = 1);
    type FaradayConstant = Real(final quantity = "FaradayConstant", final unit = "C/mol");
    type PerUnit = Real(unit = "1");
  end SIunits;
  annotation(version = "3.2.3", versionBuild = 4, versionDate = "2019-01-23", dateModified = "2020-06-04 11:00:00Z");
end Modelica;

package PL_Lib "Perpetual Library"
  extends PL_Lib.Icons.PL_icon;

  package Documentation
    extends Modelica.Icons.Information;
  end Documentation;

  package Interfaces
    extends Modelica.Icons.InterfacesPackage;

    partial model HeatExchangerBase
      outer ThermoPower.System system "System wide properties";
      replaceable package ColdFluid = Modelica.Media.Interfaces.PartialMedium annotation(choicesAllMatching = true);
      replaceable package HotFluid = Modelica.Media.Interfaces.PartialMedium annotation(choicesAllMatching = true);
      ThermoPower.Gas.FlangeA infl_1(redeclare package Medium = ColdFluid) "Cold fluid inlet";
      ThermoPower.Gas.FlangeA infl_2(redeclare package Medium = HotFluid) "Hot fluid inlet";
      ThermoPower.Gas.FlangeB outfl_1(redeclare package Medium = ColdFluid) "Cold fluid outlet";
      ThermoPower.Gas.FlangeB outfl_2(redeclare package Medium = HotFluid) "Hot fluid outlet";
    end HeatExchangerBase;

    partial model CompressorBase
      replaceable package Medium = Modelica.Media.Interfaces.PartialMedium annotation(choicesAllMatching = true);
      parameter Boolean explicitIsentropicEnthalpy = true "isentropicEnthalpy function used";
      parameter Modelica.SIunits.PerUnit eta_mech = 0.98 "mechanical efficiency";
      parameter Medium.AbsolutePressure pstart_in "inlet start pressure";
      parameter Medium.AbsolutePressure pstart_out "outlet start pressure";
      parameter Medium.Temperature Tdes_in "inlet design temperature";
      parameter Boolean allowFlowReversal = system.allowFlowReversal "= true to allow flow reversal, false restricts to design direction" annotation(Evaluate = true);
      outer ThermoPower.System system "System wide properties";
      parameter Medium.Temperature Tstart_in = Tdes_in "inlet start temperature";
      parameter Medium.Temperature Tstart_out "outlet start temperature";
      parameter Medium.MassFraction[Medium.nX] Xstart = Medium.reference_X "start gas composition";
      Medium.BaseProperties gas_in(p(start = pstart_in), T(start = Tstart_in), Xi(start = Xstart[1:Medium.nXi]));
      Medium.BaseProperties gas_iso(p(start = pstart_out), T(start = Tstart_out), Xi(start = Xstart[1:Medium.nXi]));
      Medium.SpecificEnthalpy hout_iso "Outlet isentropic enthalpy";
      Medium.SpecificEnthalpy hout "Outlet enthaply";
      Medium.SpecificEntropy s_in "Inlet specific entropy";
      Medium.AbsolutePressure pout(start = pstart_out) "Outlet pressure";
      Medium.MassFlowRate w "Gas flow rate";
      Modelica.SIunits.Angle phi "shaft rotation angle";
      Modelica.SIunits.AngularVelocity omega "shaft angular velocity";
      Modelica.SIunits.Torque tau "net torque acting on the compressor";
      Modelica.SIunits.PerUnit eta "isentropic efficiency";
      Modelica.SIunits.PerUnit PR "pressure ratio";
      ThermoPower.Gas.FlangeA inlet(redeclare package Medium = Medium, m_flow(min = if allowFlowReversal then -Modelica.Constants.inf else 0));
      ThermoPower.Gas.FlangeB outlet(redeclare package Medium = Medium, m_flow(max = if allowFlowReversal then +Modelica.Constants.inf else 0));
      Modelica.Mechanics.Rotational.Interfaces.Flange_a shaft_a;
      Modelica.Mechanics.Rotational.Interfaces.Flange_b shaft_b;
    equation
      w = inlet.m_flow;
      assert(w >= 0, "The compressor model does not support flow reversal");
      inlet.m_flow + outlet.m_flow = 0 "Mass balance";
      gas_in.p = inlet.p;
      gas_in.h = inStream(inlet.h_outflow);
      gas_in.Xi = inStream(inlet.Xi_outflow);
      outlet.p = pout;
      outlet.h_outflow = hout;
      outlet.Xi_outflow = gas_in.Xi;
      inlet.h_outflow = inStream(outlet.h_outflow);
      inlet.Xi_outflow = inStream(outlet.Xi_outflow);
      gas_iso.Xi = gas_in.Xi;
      if explicitIsentropicEnthalpy then
        hout_iso = Medium.isentropicEnthalpy(outlet.p, gas_in.state) "Approximated isentropic enthalpy";
        hout - gas_in.h = 1 / eta * (hout_iso - gas_in.h);
        s_in = 0;
        gas_iso.p = 1e5;
        gas_iso.T = 300;
      else
        gas_iso.p = pout;
        s_in = Medium.specificEntropy(gas_in.state);
        s_in = Medium.specificEntropy(gas_iso.state);
        hout - gas_in.h = 1 / eta * (gas_iso.h - gas_in.h);
        hout_iso = 0;
      end if;
      w * (hout - gas_in.h) = tau * omega * eta_mech "Energy balance";
      PR = pout / gas_in.p "Pressure ratio";
      shaft_a.phi = phi;
      shaft_b.phi = phi;
      shaft_a.tau + shaft_b.tau = tau;
      der(phi) = omega;
    end CompressorBase;
  end Interfaces;

  package Templates
    extends Modelica.Icons.BasesPackage;

    package PACK "Templates of the Passenger Air Conditioning Unit (PACK)"
      extends Modelica.Icons.BasesPackage;

      partial model PACK_simpleTemplate
        parameter Modelica.SIunits.MassFlowRate whex_cold "nominal (and initial) mass flow rate";
        parameter Modelica.SIunits.MassFlowRate whex_hot "nominal (and initial) mass flow rate";
        replaceable package HotFluid = Modelica.Media.Interfaces.PartialMedium annotation(choicesAllMatching = true);
        replaceable package ColdFluid = Modelica.Media.Interfaces.PartialMedium annotation(choicesAllMatching = true);
        replaceable Interfaces.CompressorBase compressor(redeclare package Medium = HotFluid, pstart_in = 100000, pstart_out = 100000, Tdes_in = 573.15, Tstart_out = 573.15) annotation(choicesAllMatching = true);
        replaceable Interfaces.HeatExchangerBase PHX(redeclare package ColdFluid = ColdFluid, redeclare package HotFluid = HotFluid) annotation(choicesAllMatching = true);
        ThermoPower.Gas.SourcePressure sourceP_BAin(redeclare package Medium = HotFluid, use_in_T = true, use_in_p0 = true);
        ThermoPower.Gas.SinkPressure sinkP_RA_PHXout(redeclare package Medium = ColdFluid, use_in_p0 = true);
        ThermoPower.Gas.SourceMassFlow sourceMassFlow_RA_PHXin(redeclare package Medium = ColdFluid, w0 = whex_cold, use_in_w0 = false, use_in_T = true);
        ThermoPower.Gas.SinkPressure sinkP_PACKout(redeclare package Medium = HotFluid);
      protected
        parameter Modelica.SIunits.Inertia J_shaft = 200;
        parameter Modelica.SIunits.AngularVelocity w0 = 523.3;
      equation
        connect(sourceMassFlow_RA_PHXin.flange, PHX.infl_1);
        connect(sourceP_BAin.flange, PHX.infl_2);
        connect(PHX.outfl_2, compressor.inlet);
        connect(PHX.outfl_1, sinkP_RA_PHXout.flange);
        connect(compressor.outlet, sinkP_PACKout.flange);
      end PACK_simpleTemplate;
    end PACK;
  end Templates;

  package Components
    extends Modelica.Icons.VariantsPackage;

    model HX_1DCoFlow_mass
      extends PL_Lib.Interfaces.HeatExchangerBase;
      extends PL_Lib.Icons.HeatExchanger_icon;
      import SI = Modelica.SIunits;
      parameter Integer Nnodes = 10 "number of Nodes";
      parameter Integer Nt = 20 "Number of tubes in parallel";
      parameter Real Cfhex = 0.005 "friction coefficient";
      parameter Modelica.SIunits.Length Lhex = 10 "total length";
      parameter Modelica.SIunits.Diameter Dihex = 0.02 "internal diameter";
      parameter Modelica.SIunits.Radius rhex = Dihex / 2 "internal radius";
      parameter Modelica.SIunits.Length omegahex = Modelica.Constants.pi * Dihex "internal perimeter";
      parameter Modelica.SIunits.Area Ahex = Modelica.Constants.pi * rhex ^ 2 "internal cross section";
      parameter ThermoPower.Choices.Init.Options initOpt = ThermoPower.Choices.Init.Options.steadyState "Initialisation option";
      parameter SI.AbsolutePressure pstart_c = 1e5 "Pressure start value - cold side";
      parameter SI.Temperature Tstartbar_c = 300 "Avarage temperature start value - cold side";
      parameter SI.Temperature Tstartin_c = HX_coldSide.Tstartbar "Inlet temperature start value - cold side";
      parameter SI.Temperature Tstartout_c = HX_coldSide.Tstartbar "Outlet temperature start value - cold side";
      parameter SI.AbsolutePressure pstart_h = 1e5 "Pressure start value - hot side";
      parameter SI.Temperature Tstartbar_h = 300 "Avarage temperature start value - hot side";
      parameter SI.Temperature Tstartin_h = HX_hotSide.Tstartbar "Inlet temperature start value - hot side";
      parameter SI.Temperature Tstartout_h = HX_hotSide.Tstartbar "Outlet temperature start value - hot side";
      parameter SI.Length rint = 0.01 / 2 "Internal radius (single tube)";
      parameter SI.Length rext = 0.012 / 2 "External radius (single tube)";
      parameter SI.HeatCapacity rhomcm = 7800 * 650 "Metal heat capacity per unit volume [J/m^3.K]";
      parameter SI.ThermalConductivity lambda = 20 "Thermal conductivity";
      parameter SI.MassFlowRate wnom_c = 0.1 "Nominal mass flowrate (total) - cold side";
      parameter SI.MassFlowRate wnom_h = 0.1 "Nominal mass flowrate (total) - hot side";
      parameter SI.Temperature Tstartbar_wall = 300 "Avarage temperature - wall";
      parameter SI.Temperature Tstart1_w = metalTubeFV.Tstartbar "Temperature start value - first volume - wall";
      parameter SI.Temperature TstartN_w = metalTubeFV.Tstartbar "Temperature start value - last volume - wall";
      parameter SI.Density rhohex = 1000 "Density of the material";
      parameter SI.Mass mhex = 1 "Compressor mass calculated from material properties and sizing";
      ThermoPower.Gas.Flow1DFV HX_hotSide(redeclare package Medium = HotFluid, A = Ahex, omega = omegahex, wnom = 0.1, Cfnom = Cfhex, Dhyd = Dihex, FFtype = ThermoPower.Choices.Flow1D.FFtypes.Cfnom, L = Lhex, N = Nnodes, dpnom = 1000, pstart = pstart_h, Tstartbar = Tstartbar_h, Tstartin = Tstartin_h, Tstartout = Tstartout_h, initOpt = initOpt, noInitialPressure = false);
      ThermoPower.Gas.Flow1DFV HX_coldSide(redeclare package Medium = ColdFluid, A = Ahex, omega = omegahex, wnom = 0.1, Cfnom = Cfhex, Dhyd = Dihex, FFtype = ThermoPower.Choices.Flow1D.FFtypes.Cfnom, L = Lhex, N = Nnodes, dpnom = 1000, pstart = pstart_c, Tstartbar = Tstartbar_c, Tstartin = Tstartin_c, Tstartout = Tstartout_c, initOpt = initOpt, noInitialPressure = false);
      ThermoPower.Thermal.MetalTubeFV metalTubeFV(L = Lhex, Nt = Nt, Nw = Nnodes - 1, lambda = lambda, rext = rext, rhomcm = rhomcm, rint = rint, Tstartbar = Tstartbar_wall, Tstart1 = Tstart1_w, TstartN = TstartN_w);
      ThermoPower.Thermal.HeatExchangerTopologyFV heatExchangerTopologyFV(Nw = Nnodes - 1);
    equation
      connect(infl_1, HX_coldSide.infl);
      connect(HX_coldSide.outfl, outfl_1);
      connect(infl_2, HX_hotSide.infl);
      connect(HX_hotSide.outfl, outfl_2);
      connect(HX_hotSide.wall, metalTubeFV.ext);
      connect(metalTubeFV.int, heatExchangerTopologyFV.side2);
      connect(heatExchangerTopologyFV.side1, HX_coldSide.wall);
    end HX_1DCoFlow_mass;

    model Compressor_noMaps_mass
      extends ThermoPower.Icons.Gas.Compressor;
      extends PL_Lib.Interfaces.CompressorBase;
      import ThermoPower.Choices.TurboMachinery.TableTypes;
      parameter Modelica.SIunits.AngularVelocity Ndesign = 523.3 "Design velocity";
      parameter Real[:, :] tablePhic = tablePhicC "Table for phic(N_T,beta)";
      parameter Real[:, :] tableEta = tableEtaC "Table for eta(N_T,beta)";
      parameter Real[:, :] tablePR = tablePRC "Table for eta(N_T,beta)";
      parameter String fileName = "noName" "File where matrix is stored";
      parameter TableTypes Table = TableTypes.matrix "Selection of the way of definition of table matrix";
      parameter Real eta_set = 0.95;
      parameter Real PR_set = 2.5;
      parameter Modelica.SIunits.Mass mass = 20 "Compressor mass";
      Modelica.Blocks.Tables.CombiTable2D Eta(tableOnFile = if Table == TableTypes.matrix then false else true, table = tableEta, tableName = if Table == TableTypes.matrix then "NoName" else "tabEta", fileName = if Table == TableTypes.matrix then "NoName" else fileName, smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative);
      Modelica.Blocks.Tables.CombiTable2D PressRatio(tableOnFile = if Table == TableTypes.matrix then false else true, table = tablePR, tableName = if Table == TableTypes.matrix then "NoName" else "tabPR", fileName = if Table == TableTypes.matrix then "NoName" else fileName, smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative);
      Modelica.Blocks.Tables.CombiTable2D Phic(tableOnFile = if Table == TableTypes.matrix then false else true, table = tablePhic, tableName = if Table == TableTypes.matrix then "NoName" else "tabPhic", fileName = if Table == TableTypes.matrix then "NoName" else fileName, smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative);
      Real N_T "Referred speed ";
      Real N_T_design "Referred design velocity";
      Real phic "Flow number ";
      Real beta(start = integer(size(tablePhic, 1) / 2)) "Number of beta line";
    protected
      parameter Real[6, 4] tableEtaC = [0, 95, 100, 105; 1, 82.5e-2, 81e-2, 80.5e-2; 2, 84e-2, 82.9e-2, 82e-2; 3, 83.2e-2, 82.2e-2, 81.5e-2; 4, 82.5e-2, 81.2e-2, 79e-2; 5, 79.5e-2, 78e-2, 76.5e-2];
      parameter Real[6, 4] tablePhicC = [0, 95, 100, 105; 1, 38.3e-3, 43e-3, 46.8e-3; 2, 39.3e-3, 43.8e-3, 47.9e-3; 3, 40.6e-3, 45.2e-3, 48.4e-3; 4, 41.6e-3, 46.1e-3, 48.9e-3; 5, 42.3e-3, 46.6e-3, 49.3e-3];
      parameter Real[6, 4] tablePRC = [0, 95, 100, 105; 1, 22.6, 27, 32; 2, 22, 26.6, 30.8; 3, 20.8, 25.5, 29; 4, 19, 24.3, 27.1; 5, 17, 21.5, 24.2];
    equation
      N_T_design = Ndesign / sqrt(Tdes_in) "Referred design velocity";
      N_T = 100 * omega / (sqrt(gas_in.T) * N_T_design) "Referred speed definition, as percentage of design velocity";
      phic = w * sqrt(gas_in.T) / gas_in.p "Flow number definition";
      Phic.u1 = beta;
      Phic.u2 = N_T;
      phic = Phic.y;
      Eta.u1 = beta;
      Eta.u2 = N_T;
      eta = eta_set;
      PressRatio.u1 = beta;
      PressRatio.u2 = N_T;
      PR = PR_set;
    end Compressor_noMaps_mass;
  end Components;

  package Icons
    extends Modelica.Icons.IconsPackage;

    partial package PL_icon
      extends Modelica.Icons.Package;
    end PL_icon;

    model HeatExchanger_icon end HeatExchanger_icon;
  end Icons;

  package Configurations
    extends Modelica.Icons.Package;

    model PACK_simple
      extends Templates.PACK.PACK_simpleTemplate(redeclare Components.HX_1DCoFlow_mass PHX, redeclare Components.Compressor_noMaps_mass compressor);
    end PACK_simple;
  end Configurations;

  package Experiments
    extends Modelica.Icons.ExamplesPackage;

    model SimplePack
      inner ThermoPower.System system;
      replaceable package HotFluid = Modelica.Media.Air.DryAirNasa constrainedby Modelica.Media.Interfaces.PartialMedium;
      replaceable package ColdFluid = Modelica.Media.Air.DryAirNasa constrainedby Modelica.Media.Interfaces.PartialMedium;
      extends PL_Lib.Configurations.PACK_simple(redeclare package ColdFluid = ColdFluid, redeclare package HotFluid = HotFluid);
    end SimplePack;
  end Experiments;
  annotation(version = "1.0.0");
end PL_Lib;

model SimplePack_total
  extends PL_Lib.Experiments.SimplePack;
end SimplePack_total;
