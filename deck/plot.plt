reset
###size 3.0 2.4
set term epslatex size 3.0,3.4 standalone color colortext 8

set style line 1 lc rgb 'blue' lt 1 lw 1.5 pt 7 pointsize 0.3 # --- blue
set style line 2 lc rgb 'red' lt 1 lw 1.5 pt 7 pointsize 0.3 # --- red
set style line 3 lc rgb 'green' lt 1 lw 1.5 pt 7 pointsize 0.3 # --- green
set style line 4 lc rgb 'orange' lt 1 lw 1.5 pt 7 pointsize 0.3 # --- orange
set style line 5 lc rgb 'black' lt 1 lw 1.5 pt 7 pointsize 0.3 # --- black

set style line 6 lc rgb 'blue' lt 11 lw 1.5 pt 7 pointsize 0.3 # --- blue
set style line 7 lc rgb 'red' lt 11 lw 1.5 pt 7 pointsize 0.3 # --- red
set style line 8 lc rgb 'green' lt 11 lw 1.5 pt 7 pointsize 0.3 # --- green
set style line 9 lc rgb 'orange' lt 11 lw 1.5 pt 7 pointsize 0.3 # --- orange
set style line 10 lc rgb 'black' lt 11 lw 1.5 pt 7 pointsize 0.3 # --- black

set output sprintf('%s.tex',fileout)

unset border
set border 3

set tmargin 3.0
set rmargin -1.5
set bmargin 4.0

#set style rect fc lt -1 fs solid 0.15 noborder
#set obj rect from 1.8,1.9 to 7.1,2.1 behind


list = system("ls out_deck*")

#print list

set xrange[0.8:2.5]
##set xtics (2, 3, 4, 5, 6, 7)
set xtics nomirror

##set yrange[1.5:5.1]
##set ytics (2, 3, 4, 5, 6)
set ytics nomirror


set label "\\small{$m_{\\rho\\pi}$}" at screen 0.82,0.05
##set label "\\small{$\\left | B \\right |$}" at screen 0.06,0.92
set label "\\small{$\\left | B \\right |$}" at screen 0.06,0.96
##set label sprintf("$\\mathbf{d}=[%s]$",number) at screen 0.45,0.95

p for [file in list] file u 1:2 w l notitle ##ls 10