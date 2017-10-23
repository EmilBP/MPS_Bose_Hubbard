#gnuplot file
#plot for [i=4:10] "< paste ../condensateData.txt ../Mycorrelation.txt" u 1:(abs(column(i-2)-column(i+6))) w l title ' '.i.' particles'

plot '../condensateData.txt' u 1:2 w l title '10 particles', \
     '../condensateData.txt' u 1:3 w l title '20 particles',\
     '../condensateData.txt' u 1:4 w l title '30 particles', \
     '../condensateData.txt' u 1:5 w l title '40 particles', \
     '../condensateData.txt' u 1:6 w l title '50 particles', \
     '../condensateData.txt' u 1:7 w l title '75 particles', \
     '../condensateData.txt' u 1:8 w l title '100 particles'

set key
#set logscale y
#set format y "$%.0t \\cdot 10^{%T}$"
#set format y "$10^{%T}$"
set xlabel 'U/J'
set ylabel '$\rho$'
set terminal epslatex size 12cm, 7cm color colortext standalone
set output './fig.tex'
replot
