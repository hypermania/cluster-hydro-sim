(* Local POWER_LAW source added to the moving-fluid band system.
   All symbols are dimensionless. k is frozen during a timestep. *)
ClearAll["Global`*"];
x = {rho, v, u}; dx = {drho, dv, du};
w = {rho, rho v, rho (u + epsilon v^2/2)};
birth = {rho, rho v, rho (u/2 + epsilon v^2/2)};
singleSource = -k w; binarySource = k birth;
singleJacobian = D[singleSource, {x}];
binaryJacobian = D[binarySource, {x}];
linearSource = Expand[singleSource + binarySource +
  (singleJacobian + binaryJacobian).dx];
(* Entries added to each cell's matrix/RHS; all columns belong to singles. *)
matrixSingle = -dt volume singleJacobian;
matrixBinary = -dt volume binaryJacobian;
rhsSingle = dt volume singleSource;
rhsBinary = dt volume binarySource;
checks = {
  Simplify[linearSource[[1]]] === 0,
  Simplify[linearSource[[2]]] === 0,
  Simplify[linearSource[[3]] + k (rho u + u drho + rho du)/2] === 0,
  Simplify[binaryJacobian - k {{1,0,0},{v,rho,0},
    {u/2+epsilon v^2/2,epsilon rho v,rho/2}}] === ConstantArray[0,{3,3}],
  (* Isolated donor: implicit rho sink, v and u unchanged. *)
  Simplify[(D[w,{x}]-dt singleJacobian).{-dt k rho/(1+dt k),0,0}
    -rhsSingle/volume] === {0,0,0}
};
Print["single_matrix=",InputForm[matrixSingle]];
Print["binary_matrix=",InputForm[matrixBinary]];
Print["single_rhs=",InputForm[rhsSingle]];
Print["binary_rhs=",InputForm[rhsBinary]];
Print["moving_capture_checks=",checks];
If[And@@checks,Exit[0],Exit[1]];
