****** RC-Check ******
***** Include Files *****
.OPTION POST=2
.hdl RESISTOR.va
.hdl CAPACITOR.va
.model cap simple_capacitor
.model res simple_resistor
***** Netlist *****
X1 1 2 res resistance = 2000
X2 2 0 cap capacitance = 7.9577e-8
Vs 1 0 AC 1V
***** AC Analysis *****
.AC DEC 10 1 1MEG
.print V(2)
.end