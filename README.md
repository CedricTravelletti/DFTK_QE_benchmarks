
**Important Note**: QuantumEspresso input files are very brittle. 
When adding an input file for a new test system, make sure to follow the syntax 
used in the other input files. 
In particular, make sure that sub-headers (CONTROL, ELECTRONS, ...) are left-aligned 
(no tabulation, no whitespaces). 
This can be particularly vicious, since the ASE backend used to load the system fails 
silently on such errors.
