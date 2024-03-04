(* Created in 2005 by author(s) 
Tran The Trung, 

Service d'Aeronomie,
Reduit de Verrieres - B.P.3,
Route des Gatines
91371 Verrieres le Buisson,
France

Institute of Physics,
NCST, 
46 Nguyen Van Ngoc, 
Hanoi, 
Vietnam
*)

Unprotect[I2VO1,I2VO2];
Remove[I2VO1,I2VO2];

BeginPackage["Atmosphere`I2VO`","Database`DB`"];

Unprotect[I2VO1,I2VO2];
Remove[I2VO1,I2VO2];

Off[General::spell];
Off[General::spell1];

I2VO1::usage="I2VO1[sf,trj,etau] gives the signal output in V of track 1, ODS model at Ouagadougou, for the list of solar trajectories trj (in the form {{mu1,phi1},...}), from input scattered field sf. sf is interpolation of solar angle and 2 viewing angles; etau is effective optical depth of the atmosphere from which effective direct solar field is calculated. I2VO1[sf,trj,etau,temp] gives result with temperature correction, where temp is list of temperature (in Celcius degree) with same length as trj.";

I2VO2::usage="I2VO2[sf,trj,etau] gives the signal output in V of track 2, ODS model at Ouagadougou, for the list of solar trajectories trj (in the form {{mu1,phi1},...}), from input scattered field sf. sf is interpolation of solar angle and 2 viewing angles; etau is effective optical depth of the atmosphere from which effective direct solar field is calculated. I2VO2[sf,trj,etau,temp] gives result with temperature correction, where temp is list of temperature (in Celcius degree) with same length as trj.";

I2VO1::temp="Input temperature correction has length `1` not equal to length of solar trajectory positions `2`.";
I2VO2::temp="Input temperature correction has length `1` not equal to length of solar trajectory positions `2`.";

Begin["`Private`"];

(* All in SI units, including angle in radian,
 except wavelength in nm*)

(* read in all parameters and pre-calculated results*)

(* solar spectrum outside atmosphere *)
SolSpec = ReadDB[
ToFileName[{$TopDirectory
,"AddOns","Applications","Atmosphere"},"SunSpectrum.db"]
];
SolSpec = ListInterpolation[SolSpec[[1]], SolSpec[[2]], InterpolationOrder -> 1];

(* dome transmission *)
DomeTrans = Interpolation[Import[
ToFileName[{$TopDirectory
,"AddOns","Applications","Atmosphere"},"Dome.hdf"]
,"HDF"], InterpolationOrder -> 1];

(* photodiode Hamamatsu S2386 response *)
S2386 = ReadDB[
ToFileName[{$TopDirectory
,"AddOns","Applications","Atmosphere"},"S2386.db"]
];
S2386 = ListInterpolation[S2386[[1]], S2386[[2]], InterpolationOrder -> 1];

(* blue filters *)
ug1=Transpose[Import[
ToFileName[{$TopDirectory
,"AddOns","Applications","Atmosphere"},"UG1.txt"]
,"Table"]];
bg12=Transpose[Import[
ToFileName[{$TopDirectory
,"AddOns","Applications","Atmosphere"},"BG12.txt"]
,"Table"]];
ug1 = ListInterpolation[ug1[[-1]], {ug1[[1]]},
InterpolationOrder -> 1];
bg12 = ListInterpolation[bg12[[-1]], {bg12[[1]]},
InterpolationOrder -> 1];

(* red filters*)
RedFilter = ReadDB[
ToFileName[{$TopDirectory
,"AddOns","Applications","Atmosphere"},"RedFilter.db"]
];
RedFilter = ListInterpolation[0.01RedFilter[[1]], RedFilter[[2]],
InterpolationOrder -> 1];
GrayFilter = ReadDB[
ToFileName[{$TopDirectory
,"AddOns","Applications","Atmosphere"},"GrayFilter.db"]
];
GrayFilter = ListInterpolation[0.01GrayFilter[[1]], GrayFilter[[2]],
InterpolationOrder -> 1];

(* Optical depth of Atmospheric Background *)
BackTau = Interpolation[Import[
ToFileName[{$TopDirectory
,"AddOns","Applications","Atmosphere"},"Tau.hdf"]
,"HDF"], InterpolationOrder -> 1];
Tau1 = BackTau[370];
Tau2 = BackTau[870];

(* ODS electronic transfer functions *)
i2v1 = Import[
ToFileName[{$TopDirectory
,"AddOns","Applications","Atmosphere"},"I2VO1.txt"]
, "Expression"];
i2v2 = Import[
ToFileName[{$TopDirectory
,"AddOns","Applications","Atmosphere"},"I2VO2.txt"]
, "Expression"];
I2V1[i_, Temp_] := i2v1[Temp, Log[10, i]];
I2V2[i_, Temp_] := i2v2[Temp, Log[10, i]];

(* Important:
phi angle must not including 2Pi because 
this angle is already accounted at 0
*)
AngleList60 = Import[
ToFileName[{$TopDirectory
,"AddOns","Applications","Atmosphere"},"AngleList60.hdf"]
,"HDF"];

(* 
field of view for Ouagadougou model 
xoay 94 °
Using angle list 60, 83
*)

FOVO1=Module[{
f =Import[
ToFileName[{$TopDirectory
,"AddOns","Applications","Atmosphere"},"FOV-O1.hdf"]
, "HDF"],end},
f=Select[f,(#[[2]]!=270)&];
f=Map[{#[[1]],If[#[[2]] > 266, #[[2]] - 266, 94 + #[[2]]],#[[3]]}&,f]; (* 266 = 360 - 94 *)
end = Select[f,(#[[2]]==360)&];
end = Map[{#[[1]],0,#[[3]]}&,end];
f=Join[f,end];
f=Map[{#[[1]]/180 Pi,#[[2]]/180 Pi,#[[3]]}&,f];
Interpolation[f,InterpolationOrder->1,PeriodicInterpolation->{False,True}]
];

FOVO2=Module[{
f =Import[
ToFileName[{$TopDirectory
,"AddOns","Applications","Atmosphere"},"FOV-O2.hdf"]
, "HDF"],end},
f=Select[f,(#[[2]]!=270)&];
f=Map[{#[[1]],If[#[[2]] > 266, #[[2]] - 266, 94 + #[[2]]],#[[3]]}&,f]; (* 266 = 360 - 94 *)
end = Select[f,(#[[2]]==360)&];
end = Map[{#[[1]],0,#[[3]]}&,end];
f=Join[f,end];
f=Map[{#[[1]]/180 Pi,#[[2]]/180 Pi,#[[3]]}&,f];
Interpolation[f,InterpolationOrder->1,PeriodicInterpolation->{False,True}]
];

(* 
Very important here!!! 
AngleAreaList60 must be replaced by 
AngleCosAreaList60
because we integrate on the plane of ODS pupil!!!!!
*)
SFOVO1=Module[{
aal = Import[
ToFileName[{$TopDirectory
,"AddOns","Applications","Atmosphere"},"AngleCosAreaList60.hdf"]
,"HDF"]},
Map[(FOVO1@@#)&, AngleList60] aal
];

SFOVO2=Module[{
aal = Import[
ToFileName[{$TopDirectory
,"AddOns","Applications","Atmosphere"},"AngleCosAreaList60.hdf"]
,"HDF"]},
Map[(FOVO2@@#)&, AngleList60] aal
];

(*
Show[Graphics[Map[{Hue[FOVO2@@#],
Point[#[[1]]{Cos[#[[2]]],Sin[#[[2]]]}]}&,
AngleList60]],AspectRatio->1];

Show[Graphics[MapThread[{Hue[1000 #1],
Point[#2[[1]]{Cos[#2[[2]]],Sin[#2[[2]]]}]}&,
{SFOVO2,AngleList60}]],AspectRatio->1];
*)

(* Useful flux *)
Off[NIntegrate::slwcon];
Off[NIntegrate::ncvb];
FluxODSO1=
NIntegrate[SolSpec[l]ug1[l]bg12[l]S2386[l]DomeTrans[l], {l, 300, 1100}]*
(1.5 10^(-3))^2;
FluxODSO2=
NIntegrate[SolSpec[l]GrayFilter[l]RedFilter[l]S2386[l]DomeTrans[l], {l, 300, 1100}]*
(1.5 10^(-3))^2;
On[NIntegrate::ncvb];
Off[NIntegrate::slwcon];

(*
I2VO1[{inten_,lat_,ze_,ph_},
mov:{{_?((0<=#<Pi)&),_?NumericQ}..},dtau_:0]:=Module[{si,di},
(* scattered current *)
si = Outer[(inten/.{lat->#1[[1]],ze->#2[[1]],ph->(#2[[2]]-#1[[2]])})&,
mov,AngleList60,1].SFOVO1;
(* direct current *)
di = Map[If[#[[1]]<N[90°],(FOVO1@@#)*
Pi Cos[#[[1]]] Exp[-(dtau+Tau1)/Cos[#[[1]]]],0]&,mov];
Map[I2V1[#, 40]&,FluxODSO1*(si+di)]
];

I2VO2[{inten_,lat_,ze_,ph_},
mov:{{_?((0<=#<Pi)&),_?NumericQ}..},dtau_:0]:=Module[{si,di},
(* scattered current *)
si = Outer[(inten/.{lat->#1[[1]],ze->#2[[1]],ph->(#2[[2]]-#1[[2]])})&,
mov,AngleList60,1].SFOVO2;
(* direct current *)
di = Map[If[#[[1]]<N[90°],(FOVO2@@#)*
Pi Cos[#[[1]]] Exp[-(dtau+Tau2)/Cos[#[[1]]]],0]&,mov];
Map[I2V2[#, 15]&,FluxODSO2*(si+di)]
];
*)

(* original *)
I2VO1[inten_,
mov:{{_?((0<=#<Pi)&),_?NumericQ}..},dtau_:Tau1]:=Module[{si,di},
(* scattered current *)
si = Outer[inten[#1[[1]],#2[[1]],#2[[2]]-#1[[2]]]&,
mov,AngleList60,1].SFOVO1;
(* direct current *)
di = Map[If[#[[1]]<N[90°],(FOVO1@@#)*
Pi Cos[#[[1]]] Exp[-dtau/Cos[#[[1]]]],0]&,mov];
Map[I2V1[#, 35]&,FluxODSO1*(si+di)]
];

(* special 
I2VO1[inten_,
mov:{{_?((0<=#<Pi)&),_?NumericQ}..},dtau_:Tau1]:=Module[{si,di},
(* scattered current *)
si = Outer[inten[#1[[1]],#2[[1]],#2[[2]]-#1[[2]]]&,
mov,AngleList60,1].SFOVO2;
(* direct current *)
di = Map[If[#[[1]]<N[90°],(FOVO2@@#)*
Pi Cos[#[[1]]] Exp[-dtau/Cos[#[[1]]]],0]&,mov];
Map[I2V1[#, 35]&,FluxODSO1*(si+di)]
];
*)

I2VO2[inten_,
mov:{{_?((0<=#<Pi)&),_?NumericQ}..},dtau_:Tau2]:=Module[{si,di},
(* scattered current *)
si = Outer[inten[#1[[1]],#2[[1]],#2[[2]]-#1[[2]]]&,
mov,AngleList60,1].SFOVO2;
(* direct current *)
di = Map[If[#[[1]]<N[90°],(FOVO2@@#)*
Pi Cos[#[[1]]] Exp[-dtau/Cos[#[[1]]]],0]&,mov];
Map[I2V2[#, 35]&,FluxODSO2*(si+di)]
];

I2VO1[inten_,
mov:{{_?((0<=#<Pi)&),_?NumericQ}..},dtau_,temp_]:=Module[{si,di},
Switch[Length[temp]==Length[mov],True,
(* scattered current *)
si = Outer[inten[#1[[1]],#2[[1]],#2[[2]]-#1[[2]]]&,
mov,AngleList60,1].SFOVO1;
(* direct current *)
di = Map[If[#[[1]]<N[90°],(FOVO1@@#)*
Pi Cos[#[[1]]] Exp[-dtau/Cos[#[[1]]]],0]&,mov];
MapThread[I2V1,{FluxODSO1*(si+di),temp}]
,_,Message[I2VO1::temp,Length[tem],Length[mov]]]
];

I2VO2[inten_,
mov:{{_?((0<=#<Pi)&),_?NumericQ}..},dtau_,temp_]:=Module[{si,di},
Switch[Length[temp]==Length[mov],True,
(* scattered current *)
si = Outer[inten[#1[[1]],#2[[1]],#2[[2]]-#1[[2]]]&,
mov,AngleList60,1].SFOVO2;
(* direct current *)
di = Map[If[#[[1]]<N[90°],(FOVO2@@#)*
Pi Cos[#[[1]]] Exp[-dtau/Cos[#[[1]]]],0]&,mov];
MapThread[I2V2,{FluxODSO2*(si+di),temp}]
,_,Message[I2VO2::temp,Length[tem],Length[mov]]]
];

End[];
Protect[I2VO1,I2VO2];
EndPackage[];

On[General::spell];
On[General::spell1];