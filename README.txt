Installation

To install QSOSIM, do the following:

In the folder where you downloaded and extracted the QSOSIM package:
(1) Use a text editor to open qsosim10.f and change the line 77. Set the variable ‘home’ to be equal to the folder you extracted the program into. Save file.
(2) Compile qsosim10.f by typing “gfortran -c qsosim10.f” to enact the changes.
(3) Run the makefile by typing “make”.
(4) Run QSOSIM10 by typing “./qsosim10”.