within Distillation;

model DistillationColumn

parameter Integer NI = 7 "number of stages" ;
parameter Integer Nc = 4 " number of compounds " ;
parameter Integer In_s = 3 " feed entering stage " ;
//V[i] or L[i] are the vapour and liquid mole flows leaving the plate respectively
Real V[NI+1], L[NI], Feed, Fvap, Fliq ; // various molar flows
Real RLiq ;// This is the reflux liquid mole flow
Real RR ; // Reflux ratio
Real x[NI,Nc](each start = 1/Nc, each min = 0, each max = 1) "mole fraction of liquid flow " ;
Real y[NI+1,Nc](each start = 1/Nc, each min = 0, each max = 1) "mole fraction of vapour flow "  ;
Real x_liq[Nc](each start = 1/Nc, each min = 0, each max = 1) " mole fraction of bottoms "  ;
Real y_vap[Nc](each start = 1/Nc, each min = 0, each max = 1) "mole fraction of distillate"  ;
Real z[Nc](each start = 1/Nc, each min = 0, each max = 1) "mole fraction of feed"  ;
Real k_c[Nc] "equlibrium constant for components" ;

matconn_sim In(Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {-248, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-250, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
matconn_sim Dist(Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {250, 316}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {250, 298}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
matconn_sim Bot(Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {250, -296}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {252, -300}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

   
equation

//input and solve specs
 
In.F = Feed ;
In.x_pc = z ;
Dist.F = Fvap ;
Dist.x_pc = y_vap ;
Bot.F = Fliq;
Bot.x_pc = x_liq ;

//DC equations: mole balance

//total mole balance
Feed = Fvap + Fliq ;
//Feed.*z[:] = Fvap.*y_vap[:] + Fliq.*x_liq[:] ; 

//condensor
V[1] = Fvap + RLiq  ;
RLiq = RR * Fvap ;
y[1,:] = y_vap[:] ;

//plate 1
V[1] + L[1] = RLiq + V[2]  ;
(V[1] - RLiq).*y[1,:] + L[1].*x[1,:] = V[2].*y[2,:] ;

// for plates exept 1 and nth 
for i in 2:(NI-1) loop
  
  if i==In_s then
    V[i] + L[i] = Feed + L[i-1] + V[i+1];
    V[i]*y[i,:] + L[i].*x[i,:] = Feed.*z[:] + L[i-1].*x[(i-1),:] + V[i+1].*y[i+1,:] ;
  else
    V[i] + L[i] = L[i-1] + V[i+1];
    V[i]*y[i,:] + L[i].*x[i,:] = L[i-1].*x[(i-1),:] + V[i+1].*y[i+1,:] ;
  end if;

end for;

//for NI th plate
V[NI] + L[NI] = L[NI-1] + V[NI+1] ;
V[NI].*y[NI,:] + L[NI].*x[NI,:] = L[NI-1].*x[NI-1,:] ;

//reboiler
//L[NI] = Fliq + V[NI+1] ;
L[NI].*x[NI,:] = Fliq.*x_liq[:] + V[NI+1].*y[NI+1,:] ;
//y[NI+1,:] = x[NI,:].*k_c[:] ;
x_liq[:] = x[NI,:] ;

//equilibrium relation
for i in 1:Nc loop
  for j in 1:NI loop
    x[j,i]*k_c[i] = y[j,i];
  end for;
end for;

for i in 1:NI loop

 sum(y[i,:]) = 1 ;

end for;

end DistillationColumn;
