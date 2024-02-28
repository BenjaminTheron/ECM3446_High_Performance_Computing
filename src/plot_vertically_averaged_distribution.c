# GNUPLOT script for plotting the vertically averaged distribution
# Based off the GNUPLOT script used for the ECM3446 Workshop 3

# Save the plot into a PNG file
set terminal png

# Specifies the name of the output file
set output "vertically_avrg_dist.png"

# Labels and ranges for the x-axis
set xlabel "x"
set xrange [0:30.0]

# Labels and ranges for the y-axis
set ylabel "u"
set yrange [0:0.5]

# Enforces an aspect ratio of 1
set size square

# Sets the aspect ratio to 1
set cbrange [0:1]

# Switchs off the key
set key off

# Plots the initial and final data
#plot "initial.dat" with line, 
plot "vert_avg_dist_final.dat" with line

# EoF
