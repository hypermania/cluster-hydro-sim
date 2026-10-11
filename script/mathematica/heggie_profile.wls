(* Fixed-mass Plummer profiles and hydrostatic homology for the Fig.4 fit. *)
ClearAll["Global`*"];
assumptions=a>0 && b>0 && mu>0;
rhoS[r_]=1/((1+mu) a^3) (1+r^2/a^2)^(-5/2);
rhoB[r_]=mu/((1+mu) (a b)^3) (1+r^2/(a b)^2)^(-5/2);
massS=Integrate[r^2 rhoS[r],{r,0,Infinity},Assumptions->assumptions];
massB=Integrate[r^2 rhoB[r],{r,0,Infinity},Assumptions->assumptions];
pressureNew[r_]=pressure[r/a]/a^4;
densityNew[r_]=density[r/a]/a^3;
massNew[r_]=mass[r/a];
hydro=FullSimplify[D[pressureNew[r],r]+densityNew[r]massNew[r]/r^2-
 (pressure'[r/a]+density[r/a]mass[r/a]/(r/a)^2)/a^5,assumptions && r>0];
checks=FullSimplify[{massS==1/(3(1+mu)),massB==mu/(3(1+mu)),
 massB/massS==mu,hydro==0},assumptions];
Print["profile checks: ",checks];
If[And@@(TrueQ /@ checks),Exit[0],Exit[1]];
