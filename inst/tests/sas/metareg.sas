/*9-3-2008 MRC-Epid JHZ*/

data gwa;
     infile cards dlm='09'x;
     format studyid $20.;
     input studyid$ beta se n rsquare no;
     study=_n_;
     t=beta/se;
     df=n-2;
     r2=t**2/(t**2+df);
     r=sqrt(r2);
     est1=se**2;
     est2=1/n;
     nn+n;
     if no<3;
cards;
"EPIC-Obesity"    	0.07381402	0.03331727	2416	0.002029175	1
"British 1958 BC"   	0.101629302	0.042536629	5631	0.001013072	1
"CoLaus"   	0.05399	0.022418	1479	0.003911563	1
"UK Blood Services 1"    	0.07497	0.04234	1456	0.002151659	1
"EPIC-Norfolk" 	0.045900507	0.0132976	15834	0.000752014	2
"MRC ELY"	0.002730949	0.03994781	1696	2.75884E-06	2
"North Finnish BC"	0.0561908	0.0265558	4830	0.000926492	2
"Oxford Biobank" 	0.1160995	0.0488065	1165	0.004841921	2
"UK Blood Services 2" 	0.04009	0.0416912	1562	0.000592381	2
"ALSPAC mothers" 	0.0278037	0.021094	6265	0.000277322	2
"Hertforshire Study"	0.069209018	0.0308808	2842	0.001765479	2
"SardiNIA" 	0.06	0.038	1412	0.001765017	2
"KORA"	0.095246	0.039896	1642	0.003463253	2
"NHS" 	0.041	0.035	2265	0.000606016	2
"PLCO/NCI" 	0.056	0.035	2238	0.001143592	2
"Dundee Controls 1"    	0.0617494	0.038827	1913	0.001321791	2
"Dundee Controls 2"    	0.0408132	0.04383	1501	0.000578103	2
"EFSOCH"    	0.17505	0.0425	1639	0.010256995	2
"DGI Controls" 	0.04172	0.0444	1503	0.000587877	2
"FUSION Controls"	-0.029	0.052	1291	0.00024123	2
"WTCCC/CAD" 	0.0117	0.038168367	1923	4.89121E-05	3
"WTCCC/HT" 	0.111	0.0386141	1974	0.004172835	3
"WTCCC/T2DM"	0.01431	0.035428571	1947	8.38718E-05	3
"YT2D-OXGN cases"   	0.0276825	0.0730028	1909	7.5396E-05	4
"Dundee Cases 1"    	0.0151262	0.0367537	1067	0.000159015	4
"Dundee Cases 2"    	0.0507757	0.0481712	617	0.001803341	4
"DGI Cases"    	-0.04738	0.04244	1543	0.000808138	4
"FUSION Cases"	-0.01	0.058	1094	2.72213E-05	4
;
proc print data=gwa (drop=rsquare study df no) noobs;
run;
PROC MIXED method = ml data=gwa;
     class study;
     model beta = / s cl;
     repeated / group = study;
     parms / parmsdata=gwa (rename=(est1=est)) eqcons = 1 to 20;
run;
PROC MIXED method = ml data=gwa;
     class study;
     model r = / s cl;
     repeated / group = study;
     parms / parmsdata=gwa (rename=(est2=est)) eqcons = 1 to 20;
run;

/*
available from
http://www.unc.edu/~ntbrewer/metareg.html
with adaptation by JHZ 12/11/2007
*/
DATA test1;
INPUT trial ln_or      est   year  year_m  cohort  design  country  mos;
CARDS;
1     0.348254219 0.121887759 2003  3     1     1     2     2
3     0.055663868 0.007594335 1999  -1    0     1     1     1
4     -0.968310623      0.000455847 2003  3     0     1     3     1
8     -0.199289901      0.001559349 2002  2     0     1     2     1
9     -0.768254653      0.078953618 2003  3     0     1     2     1
32    0.284104251 0.199718382 1998  -2    0     1     2     2
33    0.071692745 0.020341215 2001  1     0     1     2     1
39    0.435497395 0.067556427 1991  -9    0     1     1     2
40    -0.947515921      0.019046944 1996  -4    0     1     3     1
42    0.445625551 0.040364599 2000  0     0     0     1     2
46    0.249910566 0.001183328 2003  3     0     1     1     1
;
options nocenter ps=max ls=max;
PROC PRINT data=test1;
     var trial ln_or est;
run;
PROC MIXED method = ml data=test1 (keep=trial ln_or est);
class trial;
model ln_or = /s cl;
repeated / group = trial;
parms/ parmsdata=test1 eqcons = 1 to 11;
run;

PROC MIXED method = ml data=test1;
class trial country;
model ln_or = country mos /s cl;
repeated / group = trial;
parms/ parmsdata=test1 eqcons = 1 to 11;
run;
