#gnuplot file
plot './condensateData1.txt' u 1:2 w l title '4 particles', \
     './condensateData1.txt' u 1:3 w l title '5 particles',\
     './condensateData1.txt' u 1:4 w l title '6 particles', \
     './condensateData1.txt' u 1:5 w l title '7 particles', \
     './condensateData1.txt' u 1:6 w l title '8 particles', \
     './condensateData1.txt' u 1:7 w l title '9 particles', \
     './condensateData1.txt' u 1:8 w l title '10 particles', \
     './condensateData.txt' u 1:2 w l title '10 particles', \
     './condensateData.txt' u 1:3 w l title '20 particles',\
     './condensateData.txt' u 1:4 w l title '30 particles', \
     './condensateData.txt' u 1:5 w l title '40 particles', \
     './condensateData.txt' u 1:6 w l title '50 particles'

set key outside right bottom	
set xlabel 'U/J'
set ylabel '$\rho$'
set terminal epslatex size 12cm, 7cm color colortext standalone
set output './fig.tex'
replot
