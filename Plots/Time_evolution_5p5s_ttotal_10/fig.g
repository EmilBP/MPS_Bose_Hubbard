#gnuplot file

plot './TimeEvolutionFractionData5part5sitesCorr_oneUop_ttotal_10.txt' u 1:4 w l title 'U=1, J=0'

set key top left
#set logscale y
#set format y "$%.0t \\cdot 10^{%T}$"
#set format y "$10^{%T}$"
set xlabel 't'
set ylabel '$\rho$'
set terminal epslatex size 12cm, 7cm color colortext standalone
set output './fig.tex'
replot
