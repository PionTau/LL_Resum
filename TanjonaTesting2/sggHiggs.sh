#!/bin/bash

echo "The ggHiggs file is now being compiled..."

g++ -o Higgs ResumPlot.cpp complex_def.cpp CombResum.cpp FixedptResum.cpp AP.cpp Joint.cpp Borel.cpp HSum.cpp MellinFunc.cpp integration.cpp Luminosity.cpp -lgsl -lgslcblas -lm -lcomplex_bessel -lgfortran -lcuba `lhapdf-config --cflags 0 --ldflags 0`

echo "Done!"
echo "The Higgs binary has been generated."
