model SGNew

Real SG_0, SG_avx, SG_A, SG_B, SG, SG_d[1000], SG_z[10], a, b, c, d ; 
Real SG_gi_1[10], SG_gi_2[10], SG_ps[10], SG_p[10], SG_q[11] ;
Integer i ;
parameter Integer n = 10;

equation

i = 1 ;

SG = 0.8 ;
SG_0 = 0.7 ;
SG_avx = SG / SG_0 - 1 ;
SG_A = (SG_avx / 0.619) ^ 3 ;
SG_B = 3 ;

for i in 1:1000 loop

 SG_d[i] = (SG_A / SG_B * log(1 / (1 - (i-1) / 1000))) ^ (1 / SG_B) ;
 
end for;

for i in 1:10 loop

  SG_z[i] = exp(-SG_B / SG_A * SG_d[ (i-1) + (i-1) * 90 ] ^ SG_B) - exp(-SG_B / SG_A * SG_d[ (i-1) + (i) * 90] ^ SG_B) ;

end for ;

for i in 1:11 loop
  
  SG_q[i] = SG_B/SG_A * SG_d[ (i-1) + (i) * 90] ^ SG_B ;

end for;

a = 0.90042013692039768 ;
b = -0.0049102522593429 ;
c = 0.76738843430465553 ;
d = 0.98214367448668149 ;

for i in 1:10 loop

  SG_gi_1[i] = b - (b - a) * exp(-c * (SG_q[i] ^ d)) ;
  SG_gi_2[i] = b - (b - a) * exp(-c * (SG_q[i+1] ^ d)) ;
  SG_ps[i] = 1 / SG_z[i] * (SG_A / SG_B) ^ (1 / SG_B) * (SG_gi_1[i] - SG_gi_2[i]) ;   

end for;

for i in 1:10 loop

  SG_p[i] = SG_0 * ( 1 + SG_ps[i] );

end for;  


end SGNew;
