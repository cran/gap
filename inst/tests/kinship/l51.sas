/*15-1-2014 MRC-Epid JHZ*/

proc import datafile='ped51.dta' out=ped51 dbms=dta;
run;
data ped51;
     set ped51;
     array fm fid mid;
     do over fm;
        if fm=0 then fm=.;
     end;
run;
proc inbreed data=ped51 covar outcov=amatrix;
     var id fid mid;
run;
data k2;
     parm=1;
     row=_n_;
     set amatrix (drop=_type_ _panel_);
run;
options ps=max ls=max nocenter;
proc print noobs;run;
title kinship and multivariate;
proc mixed data=ped51 covtest asycov noclprint;
     class id;
     model qt = / solution covb;
     random id / type=lin(1) ldata=k2;
run;
proc glimmix data=ped51;
     class id;
     model qt = / solution covb;
     random id / type=lin(1) ldata=k2;
run;
proc mixed data=ped51 covtest asycov noclprint;
     class id;
     model bt = / solution covb;
     random id / type=lin(1) ldata=k2;
run;
proc glimmix data=ped51;
     class id;
     model bt = / dist=binomial solution covb;
     random id / type=lin(1) ldata=k2;
     output out=r resid=resid;
run;
proc print noobs;
run;
proc reg data=r;
     model resid=qt;
run;
