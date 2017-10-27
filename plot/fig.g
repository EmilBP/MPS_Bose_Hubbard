#gnuplot file

plot './TimeEvolutionFractionData5part5sitesCorr_oneUop_ttotal_10.txt' u 1:4 w l title 'FINALLY'

set key
#set logscale y
#set format y "$%.0t \\cdot 10^{%T}$"
#set format y "$10^{%T}$"
set xlabel 'U/J'
set ylabel '$\rho$'
set terminal epslatex size 12cm, 7cm color colortext standalone
set output './fig.tex'
replot
