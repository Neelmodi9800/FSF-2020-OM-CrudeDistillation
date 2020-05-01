within CrudeDistillation;

model AssayManager

parameter Real MW = 85 ;
parameter Real Tb = 400 ;
parameter Real SG = 0.85 ;
parameter Integer n = 10 ;
parameter Integer lim = 800 ;

Real MWA, MWa, MWb, MWB, MWavg;     // variables for MW
Real _MWsi[lim], MWx[lim], MWq[n + 1], MWsi[n + 1 ], MWz[n], MWasi[n], MWai[n] ;
parameter Real MW0 = 80 ;

Real SGA, SGB, SGavg;     // variables for SG
Real _SGsi[lim+1], SGx[lim+1], SGq[n + 1], SGz[n], SGasi[n], SGai[n], SGsi[n+1] ;
parameter Real SG0 = 0.7 ;
Integer SGi[n+1] ;

Real TbA, TbB, Tbavg;     // variables for Tb
Real Tbx[lim+1], _Tbsi[lim+1], Tbsi[n+1], Tbq[n + 1], Tbz[n], Tbasi[n], Tbai[n] ;
parameter Real Tb0 = 333 ;
Integer Tbi[n+1] ;

Real Pc[n], d15[n], Tc[n], AF[n], x[n] ; // varialbe for critical property calculation and mole fractions

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

for i in 0:lim loop

  SGx[i+1] = i/lim ;
  _SGsi[i+1] = (SGA/SGB * log ( 1 / (1 - (i-0.001)/lim) )) ^ (1/SGB) ;
  
end for ;  

for i in 0:n loop

  (SGsi[i+1], SGi[i+1]) = Modelica.Math.Vectors.interpolate(SGx,_SGsi,i/(n)) ;
  SGq[i+1] = SGB/SGA * SGsi[i+1] ^ SGB;

end for;

for i in 1:n loop
  
  SGz[i] = exp( -SGB/SGA * _SGsi[integer((i)*(lim)/(n+1))] ^ SGB ) - exp( -SGB/SGA * _SGsi[integer((i+1)*(lim)/(n+1))] ^ SGB ) ;
  SGasi[i] = 1/SGz[i] * (SGA/SGB) ^ (1/SGB) *( Gamma.PartialGamma( 1 + 1/SGB , SGq[i] ) - Gamma.PartialGamma( 1+1/SGB, SGq[i+1])) ; 
  SGai[i] = SG0 * ( 1 + SGasi[i] ) ;

end for; //equations for SG ends


Tbavg = ( Tb - Tb0 ) / Tb0 ; // equation starts for Tb calculation 

TbB = 1.5 ; 

TbA = (Tbavg/0.689) ^ TbB;

 Tbx[1] = 0/lim ;
  _Tbsi[1] = (TbA/TbB * log ( 1 / (1 - (0-(0-0.001)/lim)/(lim)) ) ) ^ (1/TbB) ;

for i in 1:lim loop

  Tbx[i+1] = i/lim ;
  _Tbsi[i+1] = (TbA/TbB * log ( 1 / (1 - (i-(i-0.001)*0.001)/(lim)) ) ) ^ (1/TbB) ;
  
end for ;  

for i in 0:n loop

  (Tbsi[i+1], Tbi[i+1]) = Modelica.Math.Vectors.interpolate(Tbx,_Tbsi,i/(n)) ;
  Tbq[i+1] = TbB/TbA * Tbsi[i+1] ^ TbB;

end for;

for i in 1:n loop
  
  Tbz[i] = exp( -TbB/TbA * _Tbsi[integer((i)*(lim)/(n+1))] ^ TbB ) - exp( -TbB/TbA * _Tbsi[integer((i+1)*(lim)/(n+1))] ^ TbB ) ;
  Tbasi[i] = 1/Tbz[i] * (TbA/TbB) ^ (1/TbB) *( Gamma.PartialGamma( 1 + 1/TbB , Tbq[i] ) - Gamma.PartialGamma( 1+1/TbB, Tbq[i+1])) ; 
  Tbai[i] = Tb0 * ( 1 + Tbasi[i] ) ;

end for; //equations for Tb ends

for i in 1:n loop // calcualtion for critical properties
  
  x[i] = MWz[i] ;
  d15[i] = SGai[i] ;
  Tc[i] = CriticalProp.Tc_Riazi( Tbai[i], d15[i] ) ;
  Pc[i] = CriticalProp.Pc_RiaziDaubert( Tbai[i], d15[i] ) ;
  AF[i] = CriticalProp.AcentricFactor_LeeKesler( Tc[i], Pc[i], Tbai[i] ) ;
  
end for;

end AssayManager;
