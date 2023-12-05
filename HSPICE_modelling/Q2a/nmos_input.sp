******** Basic nmos_test ********
***** Include Files *****
.OPTION POST=2 INGOLD=2
.hdl NMOS.va $ change the names according to your file_names
.model mos simple_NMOS $ change module name accordingly
***** NETLIST *****
X1 drain gate 0 0 mos
Vdd drain 0 DC 2.0
Vgg gate 0 DC 2.0
***** DC Analysis *****
.dc Vgg 0 2 0.02 $ to plot Input Characteristics
.print -I(Vdd) $ to get data in the .lis file
.end