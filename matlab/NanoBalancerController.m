

P = 1.75;
D = 0.15;
I=0.1;
N = 50;
s = tf('s');
Controller = P+(I/s)+D*(N/(1+(N/s)));

