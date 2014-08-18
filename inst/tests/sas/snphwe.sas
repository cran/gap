options nocenter nodate pageno=1 linesize=80 pagesize=40;

proc proto package = work.mathfun label = "SNPHWE functions";
     double SNPHWE(int obs_hets, int obs_hom1, int obs_hom2);
     link '/home/jhz22/mrc/software/SAS/snphwe.so';
run;
proc fcmp inlib=work outlib=work.mathfun.trial;
     function HWE(b,a,c);
              pHWE=SNPHWE(b,a,c);
              return(pHWE);
     endsub;
     pHWE = SNPHWE(2,1,3);
     put pHWE=;
     pHWE = SNPHWE(20,10,30);
     put pHWE=;
     pHWE = SNPHWE(200,100,300);
     put pHWE=;
run; 
options cmplib=(work work.mathfun);
data abc;
     input a b c;
     pHWE=HWE(b,a,c);
     datalines;
     1 2 3
     10 20 30
     100 200 300
     6615 9953 3774
run;
proc print;
     format pHWE 20.15;
run;

/* it is cumbersome to use "call module" from PROC IML */
