%macro mtdt(data,n);
data _bt_;
  set &data;
  array x {&n} x1-x&n;
  array allele {&n} y1-y&n;
  do i=1 to &n; allele{i}=0; end;
  y=1;
  do i=1 to &n;
     allele{_n_}=1;
     allele{i}=-1;
     count=x{i};
     if _n_ ne i then output;
     allele{i}=0;
  end;
  keep y count y1-y&n;
run;
/*Bradly-Terry model*/
proc logistic data=_bt_;
  freq count;
  model y=y1-y&n / noint;
  output out=out p=p;
  run;
/*Bowker's test of symmetry*/   
data b;
  array x x1-x&n;
  do i=1 to &n;
     set &data;
     do j=1 to &n; w=x[j]; output; end; 
  end;
  drop x1-x&n;
run;
proc freq;
   weight w;   
   table i*j / agree noprint;
   run;
%mend;
data a;
  input x1-x12;
cards;
0 0  0  2  0 0  0  0  0  0  0  0
0 0  1  3  0 0  0  2  3  0  0  0
2 3 26 35  7 0  2 10 11  3  4  1
2 3 22 26  6 2  4  4 10  2  2  0
0 1  7 10  2 0  0  2  2  1  1  0
0 0  1  4  0 1  0  1  0  0  0  0
0 2  5  4  1 1  0  0  0  2  0  0
0 0  2  6  1 0  2  0  2  0  0  0
0 3  6 19  6 0  0  2  5  3  0  0
0 0  3  1  1 0  0  0  1  0  0  0
0 0  0  2  0 0  0  0  0  0  0  0
0 0  1  0  0 0  0  0  0  0  0  0
;
%mtdt(a,12);
