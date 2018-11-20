
GenerateDdbDtg[data_] := Module[{dbLength, xs, ys, timeStamps,
   fordb, dbtk, rndStp, xdb, ydb, fortg, xtg, ytg},
  {xs, ys} = data;
  timeStamps = Range[Length[xs]];
  dbLength = IntegerPart[Length[timeStamps]*0.5] + 1;
  dbtk = Min[dbLength, 100];
  rndStp = RandomSample[timeStamps];
  fordb = Take[rndStp, dbtk];
  xdb = xs[[fordb]];
  ydb = ys[[fordb]];
  fortg = Drop[rndStp, dbtk];
  xtg = xs[[fortg]];
  ytg = ys[[fortg]];
  {xdb, ydb, xtg, ytg}]

Smap <- function(xdb, ydb, xtg, ytg, theta) := 
 Module[{yStar, ck, weights, as, bs, cs, sxdb, norms, ds},
  {yStar, ck} = Transpose[(
       sxdb = Transpose[Transpose[xdb] - INPUT_];
       norms = Norm /@ sxdb;
       ds = Mean[norms];
       weights = N[Exp[-theta*(norms/ds)]];
       as = weights*xdb;
       bs = weights*ydb;
       cs = PseudoInverse[as].bs;
       {cs.INPUT_, cs}) & /@ xtg];
  {Total[(yStar - ytg)^2]/Length[ytg],
   Mean[ck](*ctilde*)}]

(*Smap result with the best theta*)

Bestsmap[iactive_, i_, data_] := 
 Module[{activePositions, xdb, ydb, xtg, ytg,
   bestResult, mse, cTilde, theta, itest, sres, best, iscale, step, 
   neighb, newtest,
   newres},
  activePositions = Flatten[Position[iactive, 1]];
  {xdb, ydb, xtg, ytg} = {data[[1, All, activePositions]], 
    data[[2, All, i]],
    data[[3, All, activePositions]], data[[4, All, i]]};
  (*minimization of theta*)
  itest = Range[0, 8];
  sres = Sort[
    Table[{Smap[xdb, ydb, xtg, ytg, theta], theta}, {theta, itest}]];
  best = sres[[1]];
  iscale = 1; step = 0;
  While[step < 
     5 && (step != 2 || best[[2]] >= 0.5) && (step != 1 || 
      best[[2]] >= 0.5),
   (*0.5以下の場合は無視する設定.これが特にノイズが有る場合などに,
   影響しないことは確認するべき*)
   step++;
   neighb = Select[sres, (Abs[INPUT_[[2]] - best[[2]]] == iscale &)];
   newtest = (best[[2]] + neighb[[All, 2]])/2;
   newres = 
    Sort[Table[{Smap[xdb, ydb, xtg, ytg, theta], theta}, {theta, 
       newtest}]];
   sres = Sort[Join[{best}, neighb, newres]];
   best = sres[[1]];
   iscale = iscale/2
   ];
  bestResult = best;
  (*MSE and Overscript[c, ~]i (ci tilde)*)
  
  theta = N[bestResult[[2]]];
  mse = bestResult[[1, 1]];
  cTilde = bestResult[[1, 2]];
  {mse,
   ReplacePart[iactive, 
    MapThread[(INPUT_1 -> INPUT_2) &, {activePositions, cTilde}]](*cij^**), 
   theta}
  ]

SSM[timeseries_, i_, qc_, baggingIteration_] := (
  baggingResult = ParallelTable[
    numberofSpecies = Length[timeseries[[1]]];
    ddbDtg = GenerateDdbDtg[GenerateXY[timeseries]];
    iactive = 
     Append[ReplacePart[Table[0, {numberofSpecies}], i -> 1], 1];
    initialResult = Bestsmap[iactive, i, ddbDtg];
    msePrev = initialResult[[1]];
    ciStar = initialResult[[2]];
    thetaBest = initialResult[[3]];
    links = 1; q = 99;
    While[links < numberofSpecies - 1 && q > qc,
     newIactives = 
      ReplacePart[iactive, INPUT_ -> 1] & /@ 
       Flatten[Position[iactive, 0]];
     results = (Bestsmap[INPUT_, i, ddbDtg]) & /@ newIactives;
     mses = results[[All, 1]];
     mseBest = Min[mses];
     q = (1 - mseBest/msePrev);
     hBest = Position[mses, mseBest][[1, 1]];
     If[q > qc, iactive = newIactives[[hBest]]; msePrev = mseBest; 
      ciStar = results[[hBest, 2]]; thetaBest = results[[hBest, 3]]];
     links++
     ];
    {ciStar, thetaBest},
    {baggingIteration}, DistributedContexts -> Automatic, 
    Method -> "FinestGrained"];
  {Median[baggingResult[[All, 2]]],
   Drop[Median[baggingResult[[All, 1]]], -1]}
  )