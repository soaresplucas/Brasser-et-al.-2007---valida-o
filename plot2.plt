
reset
set angles degrees
set terminal pdfcairo enhanced color font "Times-Roman, 30" size 40cm,30cm
set encoding iso_8859_1
#marcação do periélio 
q1=30
q2=35
q3=40
f(x)=1-q1/x
g(x)=1-q2/x
h(x)=1-q3/x
N=3505

set colors classic
unset log
unset grid
unset key  
set key font "Times-Roman,40"
set xtics font "Times-Roman,40"
set ytics font "Times-Roman,40"
set xlabel font "Times-Roman,40"
set ylabel font "Times-Roman,40"
set angles degrees
set out 'Grafico - com_brasser_2007.pdf'
#set title "Evolução planetária"
unset label 1
unset label 2
unset label 3 
unset label 4
set multiplot layout 2,2 rowsfirst
set xlabel "t (yr)"
set ylabel " Semimajor axis (AU) "
set xrange[0:2000]
set yrange[0.94:1]
#set label 2 'J' at 1,0.6 tc ls 1

unset grid
unset key
unset log 
plot "part-1.out" u ($1):2 w l lc rgb "blue" t "e=0" , \
"part-2.out" u ($1):2 w l  lc rgb "red" t "e=0.2", \
#"part-2.out" u ($1):2 w l  lc rgb "green" t "e=0.01", \

#"part-4.out" u ($1):2 w l  lc rgb "cyan" t "e=0.3", \
#"part-5.out" u ($1):2 w l  lc rgb "pink" t "e=0.8", \
#"part-6.out" u ($1):2 w l  lc rgb "goldenrod" t "e=0.02", \
#"part-7.out" u ($1):2 w l  lc rgb "gray" t "e=0.04", \

set ylabel " Eccentricity"
set log y
set yrange[0.001:0.2]
set xrange[0:2000]
plot "part-1.out" u ($1):3 w l lc rgb "blue" t "e=0", "part-2.out" u ($1):3 w l  lc rgb "red" t "e=0.2", \

unset log
set autoscale 
set ylabel " Inclination (deg) "
unset log
set autoscale 
set log y
unset format y
unset format x
set xrange[0:2000]
set yrange[:10]
plot "part-1.out" u ($1):($4*180/pi) w l lc rgb "blue" t "e=0" , \
"part-2.out" u ($1):($4*180/pi) w l  lc rgb "red" t "e=0.2"

set key below


set ylabel "q(AU)"
unset log
set yrange[0.78:1]
set xrange[0:2000]
plot "part-1.out" u ($1):($2*(1-$3)) w l lc rgb "blue" t "e=0", "part-2.out" u ($1):($2*(1-$3)) w l  lc rgb "red" t "e=0.2", \

unset multiplot
reset


set out
set term wxt

#"corpo1.out" u ($2*cos($1)-$3*sin($1)):($2*sin($1)+$3*cos($1)) w l t "Orbita Satelite", \
#"corpo2.out" u ($2*cos($1)-$3*sin($1)):($2*sin($1)+$3*cos($1)) w l  notitle
#"corpo2.out" u ($2):($3) w l  notitle
