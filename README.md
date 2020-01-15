# HKED
A Hyper-K Event Display for WCSim

To run this code you will need to first set the WCSIMDIR environment variable. 

llib.C loads the WCSim ROOT libraries and HKED runs the event display.

You can then run the code by typing:
root -l -x llib.C 'HKED.C("wcsim.root")'
where wcsim.root is the name of your wcsim output root file. 
Alternatively just type:
./run wcsim.root
This program requires you to be running the version of WCSim with the OD
