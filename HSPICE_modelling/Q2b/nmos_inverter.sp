******** Basic NMOS Inverter ********
***** Include Files *****
.OPTION POST=2 INGOLD=2
.hdl NMOS.va
.hdl RESISTOR.va
.model mos simple_NMOS
.model res simple_resistor
***** NETLIST *****
X1 drain gate 0 0 mos
X2 supply drain res  resistance = 100000
Vdd supply 0 DC 2.0
***** Pulse Analysis ***** 
Vpulse gate 0 pulse 0 2 2n 2n 2n 98n 200n $ LV HV td tr tf PW PER
.tran 1n 300n
.option post probe
.probe V(gate) V(drain)
.print V(gate), V(drain)
.end