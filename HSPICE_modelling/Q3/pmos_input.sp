***** Basic PMOS Test *****
***** Include Files *****
.OPTION POST=2 INGOLD=2
.hdl PMOS.va $ change the names according to your file_names
.model mos simple_PMOS $ change module name accordingly

***** NETLIST *****
X1 drain gate 0 0 mos
Vdd drain 0 DC -2
Vgg gate 0 DC -2

***** DC Analysis *****
.dc Vgg 0 -2 -0.02 $ to plot Input Characteristics
.print I(Vdd) $ to get data in the .lis file
.end