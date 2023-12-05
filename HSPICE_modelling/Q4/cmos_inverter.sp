****** Basic CMOS Inverter ******
***** Include Files *****
.OPTION POST=2 INGOLD=2
.hdl NMOS.va
.hdl PMOS.va
.hdl CAPACITOR.va
.model n_mos simple_NMOS
.model p_mos simple_PMOS
.model cap simple_capacitor
***** NETLIST *****
X1 drain gate supply supply p_mos
X2 drain gate 0 0 n_mos
X3 drain 0 cap  capacitance = 50e-15
Vdd supply 0 DC 2.0
***** Pulse Analysis *****
Vpulse gate 0 pulse 0 2 2n 2n 2n 98n 200n $ LV HV td tr tf PW PER
.tran 1n 300n
.option post probe
.probe V(gate) V(drain)
.print V(gate), V(drain)
.end