(* Heggie & Aarseth 1992, Eqs. 3, 10, 15, 23, 24.
   C is denoted conductivity; every numerical solver coefficient is dimensionless. *)
ClearAll["Global`*"];
assume = G>0 && rho0>0 && r0>0 && ms>0 && mi>0 && logLambda>0;
mass0=4 Pi rho0 r0^3; energy0=G mass0/r0;
time0=r0/(Sqrt[energy0] 3 ms logLambda);
pref=FullSimplify[G^2 (ms mass0) rho0 time0/energy0^(3/2),assume];
c2=FullSimplify[24 Pi G (mi mass0) conductivity rho0 r0 Sqrt[energy0]
  time0 logLambda Sqrt[2/3]/(mass0 energy0),assume];
c1=FullSimplify[4 Sqrt[3 Pi] pref logLambda,assume];
bs=FullSimplify[(3807/1000) (3/4) pref,assume];
bb=FullSimplify[(10/4) Sqrt[3] pref,assume];
single=FullSimplify[(3807/1000) (3/2) pref logLambda fstar,assume];
relaxation=FullSimplify[(34/100) (2/3)^(3/2)/pref/logLambda,assume];
checks={pref==1/(12 Pi logLambda),c2==2 Sqrt[2/3] conductivity mi/ms,
 c1==1/Sqrt[3 Pi],bs==(3807/1000)/(16 Pi logLambda),
 bb==5 Sqrt[3]/(24 Pi logLambda),single==(3807/1000) fstar/(8 Pi),
 D[b/Sqrt[uf+uh],uf]==-b/(2(uf+uh)^(3/2)),
 D[b/Sqrt[2 uf],uf]==-b/(2 Sqrt[2] uf^(3/2))};
Print["prefactor=",InputForm[pref]," c2=",InputForm[c2]," c1=",InputForm[c1]];
Print["bs_total=",InputForm[bs]," bb=",InputForm[bb]," single=",InputForm[single]];
Print["self_relaxation_time_prefactor=",InputForm[relaxation]];
Print["heggie_unit_checks=",checks];
If[And@@(TrueQ[Simplify[#,assume && uf>0 && uh>0]]& /@checks),Exit[0],Exit[1]];
