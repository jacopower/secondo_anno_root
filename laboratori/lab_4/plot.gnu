set terminal x11

Ragg = 32.59
R = 148
L = 0.01196
C = 1.703E-7
Rtot = R + 50. + Ragg

f(x) = 5.0 * R / sqrt(Rtot * Rtot + (2 * pi * x * L - 1 / (2 * pi * x * C)) * (2 * pi * x * L - 1 / (2 * pi * x * C)))

set title "Titolo"
set xlabel "Independent Variable (no units)"
set ylabel "Dependent Variable (no units)"
set grid


plot [-1:4] f(x) title "titleDDDDF"