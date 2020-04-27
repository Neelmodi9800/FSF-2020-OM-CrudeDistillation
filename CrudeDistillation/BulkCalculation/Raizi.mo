within CrudeDistillation.BulkCalculation;

model Raizi

extends input_values;

Real MWA, MWa, MWb, MWB, MWavg;     // variables for MW
Real _MWsi[lim], MWx[lim], MWq[n + 1], MWsi[n + 1 ], MWz[n], MWasi[n], MWai[n] ;
parameter Real MW0 = 80 ;

Real SGA, SGa, SGb, SGB, SGavg;     // variables for SG
Real _SGsi[lim], SGx[lim], SGq[n + 1], SGsi[n + 1 ], SGz[n], SGasi[n], SGai[n] ;
parameter Real SG0 = 0.7 ;

Real TbA, Tba, Tbb, TbB, Tbavg;     // variables for Tb
Real _Tbsi[lim], Tbx[lim], Tbq[n + 1], Tbsi[n + 1 ], Tbz[n], Tbasi[n], Tbai[n] ;
parameter Real Tb0 = 333 ;

Integer i ;

equation

MWavg = ( MW - MW0 ) / MW0 ; // equation starts for MW calculation 

MWB = 1 ; 

MWA = (MWavg) ^ MWB;
i = 1;

for i in 1:lim loop

  MWx[i] = i/( lim + 1) ;
  _MWsi[i] = (MWA/MWB) * log ( 1 / (1 - (i)/( lim + 1)) ) ^ (1/MWB) ;
  
end for ;  

MWa = 1/n * ( MWx[lim] - MWx[1] ); // converting x as a linear function of i in 1:n
MWb = MWx[1] - MWa;

for i in 1:( n + 1) loop

  MWsi[i] = (MWA/MWB * log ( 1/ ( 1 - MWa*i - MWb ) )) ^ (1/MWB) ;
  MWq[i] = MWB/MWA * (( MWA/MWB * log ( 1/ ( 1 - MWa*i -MWb ) )) ^ (1/MWB) ) ^ MWB;

end for;

for i in 1:n loop
  
  MWz[i] = exp( -MWB/MWA * MWsi[i] ^ MWB ) - exp( -MWB/MWA * MWsi[i+1] ^ MWB ) ;
  MWasi[i] = 1/MWz[i] * (MWA/MWB) ^ (1/MWB) *( Gamma.PartialGamma( 1 + 1/MWB , MWq[i] ) - Gamma.PartialGamma( 1+1/MWB, MWq[i+1])) ;
  MWai[i] = MW0 * ( 1 + MWasi[i] );
  
end for; // MW calculation ends

SGavg = ( SG - SG0 ) / SG0 ; // equation starts for SG calculation 

SGB = 3 ; 

SGA = (SGavg/0.619) ^ SGB;

for i in 1:lim loop

  SGx[i] = i/( lim + 1) ;
  _SGsi[i] = (SGA/SGB) * log ( 1 / (1 - (i)/( lim + 1)) ) ^ (1/SGB) ;
  
end for ;  


SGa = 1/n * ( SGx[lim] - SGx[1] ); // converting x as a linear function of i in 1:n
SGb = SGx[1] - SGa;

for i in 1:( n + 1) loop

  SGsi[i] = (SGA/SGB * log ( 1/ ( 1 - SGa*i-SGb ) )) ^ (1/SGB) ;
  SGq[i] = SGB/SGA * (( SGA/SGB * log ( 1/ ( 1 - SGa * i - SGb ) )) ^ (1/SGB) ) ^ SGB;

end for;

for i in 1:n loop
  
  SGz[i] = exp( -SGB/SGA * SGsi[i] ^ SGB ) - exp( -SGB/SGA * SGsi[i+1] ^ SGB ) ;
  SGasi[i] = 1/SGz[i] * (SGA/SGB) ^ (1/SGB) *( Gamma.PartialGamma( 1 + 1/SGB , SGq[i] ) - Gamma.PartialGamma( 1+1/SGB, SGq[i+1])) ; 
  SGai[i] = SG0 * ( 1 + SGasi[i] ) ;

end for; //equations for SG ends

Tbavg = ( Tb - Tb0 ) / Tb0 ; // equation starts for Tb calculation 

TbB = 1.5 ; 

TbA = (Tbavg/0.689) ^ TbB;

for i in 1:lim loop

  Tbx[i] = i/( lim + 1) ;
  _Tbsi[i] = (TbA/TbB) * log ( 1 / (1 - (i)/( lim + 1)) ) ^ (1/TbB) ;
  
end for ;  


Tba = 1/n * ( Tbx[lim] - Tbx[1] ); // converting x as a linear function of i in 1:n
Tbb = Tbx[1] - Tba;

for i in 1:( n + 1) loop

  Tbsi[i] = (TbA/TbB * log ( 1/ ( 1 - Tba*i - Tbb ) )) ^ (1/TbB) ;
  Tbq[i] = TbB/TbA * (( TbA/TbB * log ( 1/ ( 1 - Tba*i - Tbb ) )) ^ (1/TbB) ) ^ TbB;

end for;

for i in 1:n loop
  
  Tbz[i] = exp( -TbB/TbA * Tbsi[i] ^ TbB ) - exp( -TbB/TbA * Tbsi[i+1] ^ TbB ) ;
  Tbasi[i] = 1/Tbz[i] * (TbA/TbB) ^ (1/TbB) *( Gamma.PartialGamma( 1 + 1/TbB , Tbq[i] ) - Gamma.PartialGamma( 1+1/TbB, Tbq[i+1])) ;
  Tbai[i] = Tb0 * ( 1 + Tbasi[i] ) ;
  
end for; //equations for Tb ends

end Raizi;
