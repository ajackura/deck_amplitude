#!/bin/bash

##mkdir -p deck_data
##mkdir -p deck_figures

mkdir -p deck_full_data
mkdir -p deck_full_figures

lam=0
##S=0

for M in 0 1
do
    for J in 0 1 2 #3 4
    do
	for S in 0 1 2
	do
	    ./a.out <<EOF > out_deck_J_${J}_M_${M}_S_${S}_lam_${lam}.txt
$J
$M
$S
$lam
EOF
	done
    done

    gnuplot -e "fileout='fplotM"$M"'" plot.plt
    
    latex fplotM${M}.tex

    dvips -o fplotM${M}.eps fplotM${M}.dvi

    ps2pdf fplotM${M}.eps fplotM${M}.pdf

    rm fplotM${M}.eps

    mv out_deck_J* deck_data
done

rm *.tex
rm *.dvi
rm *inc.eps
rm *.aux
rm *.log