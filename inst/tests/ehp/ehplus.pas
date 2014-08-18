{$R+,N+,E+}
{$M 65500,0,360}

PROGRAM eh(input,output);
 {This program is a utility program which estimates gene frequencies
 under the assumption of allelic independence and haplotype frequencies
 allowing for allelic associations. It also computes the Log likelihood,
 chi-square and the number of degrees of freedom.

 Updated Dec, 1993
 Added Case-Control Study
 EHIN.DAT       Required input file
 }

LABEL 999,plus;

CONST
 maxalle=30; {maximum number of alleles}
 maxloci=10; {maximum number of loci}
 maxgeno=1800;{maximum number of genotypes, ie. product over numbers of}
              {genotypes at all loci, L1*(L1+1)/2 * L2*(L2+1)/2 * ...}
 maxhap =500; {L1 * L2 * ..., where Li = number of alleles at i-th locus}
 v='1.14';  {20 Oct 97}
 bs=chr(8);
 verysmall=0.1e-8;
{/**/}
 maxposcom=400;
 maxbinsize=12;
 prime=281;
 ehplus=true;

TYPE
 real=single;
 arr=ARRAY[1..maxloci,1..maxalle] OF real;
 larr=ARRAY[1..maxloci] OF integer;
 linetype=ARRAY[1..maxhap] OF real;
 lineint=ARRAY[1..maxposcom] OF real;
{/**/}
{
 hashtable=record
   n:integer;
   obs:array[0..maxbinsize]of integer;
   id:array[0..maxbinsize]of integer;
 end;
}
VAR
 p:arr; {stores allele frequencies of each locus}
 locihn,locigeno,locik,loci,caseloci,locip,lociq:larr;  {pointer to each locus}
 i,totalloci,cnloci,haplnum,genonum:integer;
 maxall,loop,loop1:integer;
 controlf,casef,wf,tmp:text;
 totalg,oldloglike,loglike,iniloglike,diff,indloglike,oldlog:real;
 indivN,caseN,contrlN,dg11,dg12,dg22,dg33:lineint;
 oldhaplo,newhaplo,inihaplo,indephaplo:linetype;
 casestudy,done,onemark:boolean;
 fname,fnameout:STRING;
 ch:char;
 ff,casep,controlp:array[1..3] of real;
 pp,totalall,totalind,genefreq:real;
 exitsave:pointer; {specific for Borland Pascal}
{/**/}
 id:array[1..3*maxposcom]of longint;
 idsave:array[1..maxposcom]of longint;
{
 bin:array[0..prime]of hashtable;
}
 genotype,bb:boolean;
 s2:real;
 bindex,mmid,j,j1,j2,n:integer;
 numloci,sample_size,obscom:integer;
 newtime,obsid,saveid:integer;

 PROCEDURE ehlogo;
  {Displays banner for EH program}
  {Global constant used: v;  program version, char constant of length 4}
 VAR i:integer;
 BEGIN
  FOR i:=1 TO 5 DO writeln;
  writeln('        .------------------------------------------.');
  writeln('        |                                          |');
  writeln('        |          Program  EH  version ',v,'       |');
  writeln('        |                                          |');
  writeln('        `------------------------------------------''');
  writeln;
  writeln('        Programmed by Xiaoli Xie           May 1990');
  writeln;
 END;   {vplogo}


 PROCEDURE cls;
 var ii:integer;
 BEGIN
  for ii:=1 to 20 do WRITELN;
 END;

 PROCEDURE initial(VAR p:arr);
 VAR i,j:integer;

 BEGIN
  FOR i:=1 TO totalloci DO
   FOR j:=1 TO maxall DO p[i,j]:=0;
  totalind:=0;
 END;


 PROCEDURE programsuccess;
 BEGIN
  writeln; writeln;
  writeln('Program EH completed successfully. Output in file "',fnameout,'"');
 END;


 PROCEDURE programaborted;
 BEGIN
  writeln('Program EH aborted',chr(7));
 END;


 procedure getfreq;
 var sumcase,sumcontrol:real;
  i:integer;
 begin
  writeln('Enter gene frequency of disease allele');
  readln(pp);
  writeln('Enter penetrances for each genotype in the following order');
  writeln(' Genotype   +/+   +/D    D/D');
  readln(ff[1],ff[2],ff[3]);
  sumcase:=ff[1]*sqr(1-pp)+2*ff[2]*pp*(1-pp)+ff[3]*sqr(pp);
  sumcontrol:=(1-ff[1])*sqr(1-pp)+2*(1-ff[2])*pp*(1-pp)+(1-ff[3])*sqr(pp);
  casep[1]:=(ff[1]*sqr(1-pp))/sumcase;
  casep[2]:=(2*ff[2]*pp*(1-pp))/sumcase;
  casep[3]:=(ff[3]*sqr(pp))/sumcase;
  controlp[1]:=((1-ff[1])*sqr(1-pp))/sumcontrol;
  controlp[2]:=(2*(1-ff[2])*pp*(1-pp))/sumcontrol;
  controlp[3]:=((1-ff[3])*sqr(pp))/sumcontrol;
 end;

 PROCEDURE outind(var ff:text;locia:larr;a:lineint;sp:integer);
 VAR i,j,n:integer;
 BEGIN
  n:=1;
  FOR i:=sp TO totalloci DO n:=n*locia[i];
  j:=locia[totalloci];
  FOR i:=1 TO n DO
   BEGIN
    write(ff,' ',a[i]:5:2);
    IF ((i MOD j)=0) THEN writeln(ff);
   END;
 END;


 procedure casecontrol(VAR f1,f2:text;l:larr;k,totalloci:integer;VAR line:integer);
 {This procedure reads in two sets of data, case and control, and then create
  a data file based on the information from the input data for case-control study}
 VAR i,j:integer;
     ca,co:real;

 BEGIN
  FOR i:=1 TO l[k] DO
   FOR j:=1 TO i DO
    IF k=totalloci THEN
     BEGIN
      IF NOT seekeoln(f1) THEN read(f1,co);
      if not seekeoln(f2) then read(f2,ca);
      line:=line+1;
      contrlN[line]:=co;
      casen[line]:=ca;
      dg33[line]:=co+ca;
      DG11[line]:=casep[1]*ca+controlp[1]*co;
      DG12[line]:=casep[2]*ca+controlp[2]*co;
      DG22[line]:=casep[3]*ca+controlp[3]*co;
     END
    ELSE casecontrol(f1,f2,l,k+1,totalloci,line);
  IF k=totalloci THEN
  begin readln(f1); readln(f2); end;
  end;

 procedure initdg;
 var tmp:text;
     i:integer;

 begin
    {generates temp.dat file }
   assign(tmp,'temp.dat');
   rewrite(tmp);
   write(tmp,'2 ');
   for i:=1 to totalloci do write(tmp,' ',loci[i]:2);
   writeln(tmp);
   outind(tmp,locigeno,dg11,1);
   outind(tmp,locigeno,dg12,1);
   outind(tmp,locigeno,dg22,1);
   close(tmp);
 END;

 FUNCTION getN(VAR f:text;l:larr;k,totalloci:integer;VAR p:arr;VAR line:integer):real;
 {This procedure computes alleles frequencies}
 VAR i,j:integer;
     total,n,tempn,temp:real;

 BEGIN
  tempN:=0;
  FOR i:=1 TO l[k] DO
   FOR j:=1 TO i DO
    IF k=totalloci THEN
     BEGIN
      IF NOT seekeoln(f) THEN read(f,n);
      line:=line+1;
      indivN[line]:=n;
      p[k,i]:=p[k,i]+n;
      p[k,j]:=p[k,j]+n;
      tempN:=tempN+n;
     END
    ELSE BEGIN
     temp:=getN(f,l,k+1,totalloci,p,line);
     p[k,i]:=p[k,i]+temp;
     p[k,j]:=p[k,j]+temp;
     tempN:=tempN+temp;
    END;
  IF k=totalloci THEN readln(f);
  getN:=tempN;
 END;

 PROCEDURE calcula(totalloci:integer);
 {compute alleles' frequencies}
 VAR i,j:integer;
 BEGIN
  FOR i:=1 TO totalloci DO
   FOR j:=1 TO maxall DO
    p[i,j]:=p[i,j]/totalall;
 END;

 function bsch(k,n:integer):integer;
 var j,lo,hi,mid:integer;
     found:boolean;
 BEGIN
   lo:=0;hi:=n-1;
   bsch:=-1;
   found:=false;
   While((lo<=hi) and (not found)) do
   {
   Begin
     mid:=(lo+hi) div 2;
     mmid:=mid;
     if(bin[bindex].id[mid]<k) then lo:=mid+1
     else if(bin[bindex].id[mid]>k) then hi:=mid-1
     else begin bsch:=bin[bindex].obs[mid]; found:=true; end;
   end;
   }
 END;

 FUNCTION linenum(locia,locip:larr;sp:integer):integer;
 VAR sum,j:integer;

 BEGIN
  sum:=0;
  FOR j:=sp TO totalloci DO
   IF j=totalloci THEN sum:=sum+locip[j]
   ELSE sum:=(sum+(locip[j]-1))*locia[j+1];
  {/**/}
  bindex:=sum mod prime;
  {if(ehplus and genotype) then linenum:=bsch(sum,bin[bindex].n)
  else }linenum:=sum;
 END;

 procedure normalizeH(locia:larr;var a:linetype);
 { Force fixed gene freq. }
 var
  sump,sumq,asump,asumq:real;
  i,n:integer;

 begin
    sump:=0;
    sumq:=0;
    asump:=0;
    asumq:=0;
    n:=1;
    for i:=1 to totalloci do n:=n*locia[i];
    for i:=1 to (n div 2) do
    begin
      sump:=sump+inihaplo[i];
      asump:=asump+a[i];
    end;
    for i:=(n div 2) + 1 to n do
    begin
      sumq:=sumq+inihaplo[i];
      asumq:=asumq+a[i];
    end;
    for i:=1 to (n div 2) do
    begin
       inihaplo[i]:=inihaplo[i]*(1-pp)/sump;
       a[i]:=a[i]*(1-pp)/asump;
    end;
   for i:=(n div 2)+1 to n do
    begin
       inihaplo[i]:=inihaplo[i]*pp/sumq;
       a[i]:=a[i]*pp/asumq;
    end;
  end;  {NormalizeH}


 PROCEDURE outvec(locia:larr;a:linetype);
 VAR times,i,line,n:integer;

 BEGIN
  times:=1;
  n:=1;
  FOR i:=1 TO totalloci DO
   BEGIN
    n:=n*locia[i];
    locip[i]:=1;
   END;
  writeln(wf,'There are ',n,' Possible Haplotypes of These ',totalloci,' Loci.');
  writeln(wf,'They are Listed Below, with their Estimated Frequencies:');
{  writeln(wf);}
  write(wf,'-');
  FOR i:=1 TO totalloci DO write(wf,'---------');
  if casestudy and not onemark then FOR i:=1 TO 41 DO write(wf,'-')
               else for i:=1 to 33 do write(wf,'-');
  writeln(wf);
  write(wf,'|');
  FOR i:=1 TO totalloci DO write(wf,' Allele  ');
  if casestudy and not onemark
          then writeln(wf,'|         Haplotype   Frequency         |')
          else writeln(wf,'|      Haplotype Frequency      |');
  write(wf,'|');
  FOR i:=1 TO totalloci DO write(wf,'   at    ');
  if casestudy and not onemark
    then writeln(wf,'|                                       |')
    else writeln(wf,'|                               |');
  write(wf,'|');
  if casestudy then
  begin
    write(wf,' Disease ');
    for i:=1 to totalloci-1 do
    IF totalloci>9 THEN write(wf,' Marker',i:2)
                   ELSE write(wf,' Marker',i:1,' ');
  end else
  FOR i:=1 TO totalloci DO
   IF totalloci>9 THEN write(wf,' Locus ',i:2)
                  ELSE write(wf,' Locus ',i:1,' ');
  if casestudy and not onemark
    then writeln(wf,'|  Independent   Ind-Disease   w/Asso.  |')
    else writeln(wf,'|  Independent   w/Association  |');
  write(wf,'-');
  FOR i:=1 TO totalloci DO write(wf,'---------');
  if casestudy and not onemark then FOR i:=1 TO 41 DO write(wf,'-')
               else for i:=1 to 33 do write(wf,'-');
  writeln(wf);

  REPEAT
   line:=linenum(locia,locip,1);
   FOR i:=1 TO totalloci DO
   if casestudy and (i=1)
    then if locip[i]=1 then write(wf,'+':4,'     ')
                       else write(wf,'D':4,'     ')
    else write(wf,locip[i]:4,'     ');
   if casestudy and not onemark
    then writeln(wf,'      ',inihaplo[line]:8:6,'      ',indephaplo[line]:8:6,'     ',a[line]:8:6)
    else writeln(wf,'      ',inihaplo[line]:8:6,'     ',a[line]:8:6);
{    then writeln(wf,'      ',inihaplo[line],'      ',indephaplo[line],'     ',a[line])
    else writeln(wf,'      ',inihaplo[line],'     ',a[line]);}

   locip[totalloci]:=locip[totalloci]+1;
   FOR i:=totalloci DOWNTO 2 DO
    IF locip[i]>loci[i] THEN
     BEGIN locip[i]:=1;
      locip[i-1]:=locip[i-1]+1;
     END;
   times:=times+1;
  UNTIL times>n;
  write(wf,'-');
  FOR i:=1 TO totalloci DO write(wf,'---------');
  if casestudy and not onemark then writeln(wf,'-----------------------------------------')
    else writeln(wf,'---------------------------------');
  writeln(wf,'# of Iterations = ',loop);
  writeln(wf);
  writeln(wf,'                                   #param   Ln(L)     Chi-square');
  writeln(wf,'-------------------------------------------------------------------');
  if casestudy then
  begin
    n:=0;
    for i:=2 to totalloci do n:=n + (loci[i]-1);
    writeln(wf,'H0: No Association                     ',n:2,'  ',
      iniloglike:8:2,'      0.00');
    n:=1;
    for i:=2 to totalloci do n:=n*loci[i];
    n:=n-1;
    if not onemark then  {at lease one marker and one disease gene}
    begin
      if abs(indloglike-iniloglike) < verysmall then
       writeln(wf,'H1: Markers Asso., Indep. of Disease   ',n:2,'  ',
         indloglike:8:2,'      0.00')
       else writeln(wf,'H1: Markers Asso., Indep. of Disease   ',n:2,'  ',
         indloglike:8:2,(2*(indloglike-iniloglike)):10:2);
      if abs(loglike-iniloglike) < verysmall then
       writeln(wf,'H2: Markers and Disease Associated     ',(2*n):2,'  ',
         loglike:8:2,'      0.00')
       else  writeln(wf,'H2: Markers and Disease Associated     ',(2*n):2,
         '  ',loglike:8:2,(2*(loglike-iniloglike)):10:2);
    end else
      if abs(loglike-iniloglike) < verysmall then
      writeln(wf,'H1: Markers and Disease Associated     ',(2*n):2,'  ',
        loglike:8:2,'      0.00')
      else  writeln(wf,'H1: Markers and Disease Associated     ',(2*n):2,
        '  ',loglike:8:2,(2*(loglike-iniloglike)):10:2);
  end else
  begin
    n:=0;
    for i:=1 to totalloci do n:=n + (loci[i]-1);
    writeln(wf,'H0: No Association                     ',n:2,'  ',
      iniloglike:8:2,'      0.00');
    n:=1;
    for i:=1 to totalloci do n:=n*loci[i];
    n:=n-1;
    write(wf,'H1: Allelic Associations Allowed       ',n:2,'  ',loglike:8:2);
    if abs(loglike-iniloglike) < verysmall then
    writeln(wf,'      0.00')
    else writeln(wf,(2*(loglike-iniloglike)):10:2);
   end;
   writeln(wf);
 END;

 PROCEDURE outindP(totalloci:integer);
 VAR i,j:integer;
 BEGIN
  writeln(wf,'Estimates of Gene Frequencies (Assuming Independence)');
  if casestudy then writeln(wf,'(Disease gene frequencies are user specified)');
  write(wf,'----\------------');
  FOR i:=1 TO maxall DO write(wf,'--------');
  writeln(wf);
  write(wf,'locus \ allele ');
  FOR i:=1 TO maxall DO write(wf,' ':3,i:2,' ':3);
  writeln(wf);
  write(wf,'--------\--------');
  FOR i:=1 TO maxall DO write(wf,'--------');
  writeln(wf);
  if casestudy then
  begin
    write(wf,'Disease |',' ':7);
    FOR j:=1 TO maxall DO
     IF p[1,j]<verysmall THEN write(wf,' ':8)
     ELSE write(wf,p[1,j]:8:4);
     writeln(wf);

   FOR i:=2 TO totalloci DO
   BEGIN
    write(wf,(i-1):4,'    |',' ':7);
    FOR j:=1 TO maxall DO IF p[i,j]<verysmall THEN write(wf,' ':8)
     ELSE write(wf,p[i,j]:8:4);
     writeln(wf);
   END;
  end else
  FOR i:=1 TO totalloci DO
   BEGIN
    write(wf,i:4,'    |',' ':7);
    FOR j:=1 TO maxall DO IF p[i,j]<verysmall THEN write(wf,' ':8)
     ELSE write(wf,p[i,j]:8:4);
     writeln(wf);
   END;
  write(wf,'-----------------');
  FOR i:=1 TO maxall DO write(wf,'--------');
  writeln(wf);
  writeln(wf,'# of Typed Individuals: ',round(totalind):2);
  writeln(wf);
 END;


 PROCEDURE getlocih(totalloci:integer);
 {compute  the maximum number of possible genotypes}
 VAR i:integer;
 BEGIN
  FOR i:=1 TO totalloci DO locigeno[i]:=loci[i]*(loci[i]+1) DIV 2;
 END;


 PROCEDURE gethapP(totalloci:integer);
 {compute haplotype frequencies based on allele frequencies}
 VAR line,i,n,times:integer;

 BEGIN
  n:=1;
  times:=1;
  FOR i:=1 TO totalloci DO
   BEGIN locip[i]:=1;
    n:=n*loci[i];
   END;
  REPEAT
   line:=linenum(loci,locip,1);
   newhaplo[line]:=1;
   FOR i:=1 TO totalloci DO newhaplo[line]:=newhaplo[line]*p[i,locip[i]];
   locip[totalloci]:=locip[totalloci]+1;
   FOR i:=totalloci DOWNTO 2 DO
    IF locip[i]>loci[i] THEN
     BEGIN locip[i]:=1;
      locip[i-1]:=locip[i-1]+1;
     END;
   times:=times+1;
  UNTIL times>n;
  FOR i:=1 TO n DO inihaplo[i]:=newhaplo[i];
 END;  {gethapP}


 PROCEDURE getfirst(VAR first,last:integer);
  {find the first locus which is heterozygous genotype}
 VAR j:integer;
  ini:boolean;
 BEGIN
  first:=1;
  last:=1;
  ini:=true;
  FOR j:=1 TO totalloci DO
   BEGIN
    IF locihn[j]=1 THEN
     BEGIN
      IF ini THEN IF first<=j THEN BEGIN
       first:=j;
       ini:=false;
      END;
      IF j>last THEN last:=j;
     END;
   END;
 END;

 PROCEDURE swap(VAR a,b:integer);
 VAR temp:integer;
 BEGIN
  temp:=a;
  a:=b;
  b:=temp;
 END;

 PROCEDURE newtoold;
 VAR i:integer;
 BEGIN
  FOR i:=1 TO haplnum DO
   BEGIN
    oldhaplo[i]:=newhaplo[i];
    newhaplo[i]:=0;
   END;
 END;

 PROCEDURE iniindex;
 VAR i:integer;
 BEGIN
  FOR i:=1 TO totalloci DO
   BEGIN
    locip[i]:=1;
    lociq[i]:=1;
    locik[i]:=1;
    locihn[i]:=0;
   END;
 END;

 PROCEDURE inivar(nloci:integer;locii:larr);
 VAR i:integer;
 BEGIN
  haplnum:=1;
  genonum:=1;
  FOR i:=1 TO nloci DO
   BEGIN
    haplnum:=haplnum*locii[i];
    locigeno[i]:=locii[i]*(locii[i]+1) DIV 2;
    genonum:=genonum*locigeno[i];
   END;
 END;

 FUNCTION exp2(h:integer):integer;
 VAR i,temp:integer;
 BEGIN
  i:=1;
  temp:=1;
  REPEAT
   temp:=temp*2;
   i:=i+1;
  UNTIL i>h;
  exp2:=temp;
 END;


 FUNCTION inputdataok(VAR inf,outf:text;var nloci:integer;var locii:larr):boolean;
  {This function checks the input data. The first line is the number
  of alleles at each locus. Rest of the line is the numbers of individuals
  with given genotypes.
  The total number of genotypes = sum( loci(i)*(loci(i)+1)/2) i=1,2,...}

 VAR
     data:real;
     b:boolean;

 BEGIN
  nloci:=0;
  maxall:=0;
  b:=false;
  WHILE NOT seekeof(inf) DO
   IF NOT b THEN
    BEGIN
     WHILE NOT seekeoln(inf) DO
      BEGIN
       nloci:=nloci+1;
       read(inf,locii[nloci]);
       IF locii[nloci]>maxall THEN maxall:=locii[nloci];
      END;
     readln(inf);
     inivar(nloci,locii);
     writeln(outf);
     b:=true;
     i:=0;
   END ELSE
    BEGIN
     WHILE NOT seekeoln(inf) DO
      BEGIN
       read(inf,data);
       i:=i+1;
      END;
     readln(inf);
    END;
  inputdataok:=(i=genonum);
 END;

 function dataconsist(loci1,loci2:larr;n1,n2:integer):boolean;
 var
   i:integer;
   equal:boolean;

 begin
   i:=1;
   equal:=n1=n2;
   while equal and (i<=n1) do
   begin
     equal:=loci1[i]=loci2[i];
     i:=i+1;
   end;
   dataconsist:=equal;
  end;

 function likelihood(hap:linetype;totalloci:integer;dma:boolean;dg:boolean):real;
  {compute log likelihood of the data with given haplotype frequencies}
 VAR i,j,time,first,last,h:integer;
  linek,lineq,linep:integer;
  multilog,like,like11,like12,like22,temp,temp1,sa,su:real;


  PROCEDURE multilike(i:integer;VAR logtemp:real);
   {compute the haplotype frequencies on a single cell which has more than one
   genotypes}
  VAR j,k:integer;

  BEGIN {multilike}
   IF locihn[i]=1 THEN
    BEGIN
     IF i=first THEN
      IF i=last THEN
       BEGIN
        linek:=linenum(loci,locik,1);
        lineq:=linenum(loci,lociq,1);
        logtemp:=logtemp+2*hap[linek]*hap[lineq];
       END
      ELSE multilike(i+1,logtemp)
     ELSE BEGIN
      FOR j:=1 TO 2 DO
       BEGIN IF j=2 THEN swap(locik[i],lociq[i]);
        IF i=last THEN
         BEGIN
          linek:=linenum(loci,locik,1);
          lineq:=linenum(loci,lociq,1);
          logtemp:=logtemp+2*hap[linek]*hap[lineq];
         END
        ELSE multilike(i+1,logtemp);
       END;
      swap(locik[i],lociq[i]);
     END;
   END ELSE multilike(i+1,logtemp);
  END;  {multilike}


procedure probpcell(var prob:real);
var i:integer;

begin
   prob:=0;
   FOR i:=1 TO totalloci DO locihn[i]:=0;
   FOR i:=1 TO totalloci DO IF (locik[i]<>lociq[i]) THEN locihn[i]:=1;
   h:=0;
   FOR i:=1 TO totalloci DO h:=h+locihn[i];
   {/**/}
   {
   genotype:=true;
   if dma then linep:=linenum(locigeno,locip,2)
          else linep:=linenum(locigeno,locip,1);
   genotype:=false;
   }
   linep:=saveid;

   IF h<2 THEN
    BEGIN
     linek:=linenum(loci,locik,1);
     lineq:=linenum(loci,lociq,1);
     temp:=hap[linek]*hap[lineq];
     if temp > verysmall then
     begin
       if dma then
       begin
       IF h=0 THEN prob:=prob+temp {hap[linek]*hap[lineq]}
            ELSE prob:=prob+2*temp;  {hap[linek]*hap[lineq];}
       end else
       begin
       if linep <> -1 then
       IF h=0 THEN prob:=prob+ln(temp)*indivN[linep]
            ELSE prob:=prob+ln(2*temp)*indivN[linep];
       end;
     end;
    END
   ELSE BEGIN
    getfirst(first,last);
    multilog:=0;
    multilike(first,multilog);
    if multilog > verysmall then
      if dma then prob:=prob+multilog
      else if linep<>-1 then prob:=prob+ln(multilog)*indivn[linep];  {xie}
   END;
{disease +/D}
end;

 BEGIN  {likelihood}
 {compute sa and su}
 if dma then
 begin
  iniindex;
  time:=1;newtime:=1;obsid:=1;
  sa:=0;
  su:=0;
  REPEAT
{/**/}
  if ehplus then begin
     if newtime=idsave[obsid] then begin
        saveid:=obsid; obsid:=obsid+1; end else saveid:=-1;
     end else saveid:=newtime;
 { disease marker +/+ }
  lociq[1]:=1;
  locik[1]:=1;
  probpcell(like11);

{ disease marker +/D}
  lociq[1]:=2;
  locik[1]:=1;
  probpcell(like12);

 { disease marker D/D}
  locik[1]:=2;
  lociq[1]:=2;
  probpcell(like22);

  sa:=sa+ like11*ff[1] + like12*ff[2] + like22*ff[3];
  su:=su+like11*(1-ff[1])+like12*(1-ff[2])+like22*(1-ff[3]);
  {
  locip[totalloci]:=locip[totalloci]+1;
   FOR i:=totalloci DOWNTO 2 DO  (*was 2*)
    IF locip[i]>locigeno[i] THEN
     BEGIN
      locip[i]:=1;
      locip[i-1]:=locip[i-1]+1;
     END;
  }
   locik[totalloci]:=locik[totalloci]+1;
   FOR i:=totalloci DOWNTO 2 DO  {was 2}
    IF locik[i]>lociq[i] THEN
     BEGIN
      lociq[i]:=lociq[i]+1;
      locik[i]:=1;
      IF lociq[i]>loci[i] THEN
       BEGIN
        lociq[i]:=1;
        locik[i-1]:=locik[i-1]+1;
        IF i=2 THEN {was 1}
         IF locik[i-1]>lociq[i-1] THEN
          BEGIN lociq[i-1]:=lociq[i-1]+1;
           locik[i-1]:=1;
          END;
       END;
     end;
    newtime:=newtime+1;
    time:=time+3;
  UNTIL  time > genonum;     {genonum/3 for disease marker associate}
 end;

{ compute log likelihood }
  iniindex;
  time:=1;newtime:=1;obsid:=1;
  like:=0;

  REPEAT
  if ehplus then begin
     if newtime=idsave[obsid] then begin
        saveid:=obsid;obsid:=obsid+1; end else saveid:=-1;
     end else saveid:=newtime;
 { disease marker +/+ }
 if dma then
 begin
  lociq[1]:=1;
  locik[1]:=1;
  probpcell(like11);

  lociq[1]:=2;
  locik[1]:=1;
  probpcell(like12);

  locik[1]:=2;
  lociq[1]:=2;
  probpcell(like22);

 {normalized by dividing each cell probability by sa or su}
  {/**/}
  if linep<>-1 then
  begin
  if caseN[linep] > 0.01 then   like:=like+caseN[linep] * ln ( (like11*ff[1] + like12*ff[2] + like22*ff[3] )/sa) ;
  if contrlN[linep] > 0.01 then like:=like+contrlN[linep]*ln( ( like11*(1-ff[1])+like12*(1-ff[2])+like22*(1-ff[3])) /su);
  end;

  temp:=(ff[1]*like11+ff[2]*like12+ff[3]*like22);
  temp1:=((1-ff[1])*like11+(1-ff[2])*like12+(1-ff[3])*like22);

  if linep<>-1 then
  BEGIN
  dg11[linep]:=0;
  dg12[linep]:=0;
  dg22[linep]:=0;

  if casen[linep]>0.01 then begin
   dg11[linep]:=dg11[linep]+like11*casen[linep]*ff[1]/temp;
   dg12[linep]:=dg12[linep]+like12*casen[linep]*ff[2]/temp;
   dg22[linep]:=dg22[linep]+like22*casen[linep]*ff[3]/temp;end;

  if contrln[linep]>0.01 then begin
   dg11[linep]:=dg11[linep]+like11 *contrln[linep]*(1-ff[1])/temp1;
   dg12[linep]:=dg12[linep]+like12 *contrln[linep]*(1-ff[2])/temp1;
   dg22[linep]:=dg22[linep]+like22 *contrln[linep]*(1-ff[3])/temp1;end;
  END;
 end else if saveid<>-1 then
 begin
   probpcell(like11);
   like:=like + like11;
 end;
 {
 locip[totalloci]:=locip[totalloci]+1;
   FOR i:=totalloci DOWNTO 2 DO  (* was 2*)
    IF locip[i]>locigeno[i] THEN
     BEGIN
      locip[i]:=1;
      locip[i-1]:=locip[i-1]+1;
     END;
 }
   locik[totalloci]:=locik[totalloci]+1;
   FOR i:=totalloci DOWNTO 2 DO  {was 2}
    IF locik[i]>lociq[i] THEN
     BEGIN
      lociq[i]:=lociq[i]+1;
      locik[i]:=1;
      IF lociq[i]>loci[i] THEN
       BEGIN
        lociq[i]:=1;
        locik[i-1]:=locik[i-1]+1;
        IF i=2 THEN {was 1}
         IF locik[i-1]>lociq[i-1] THEN
          BEGIN lociq[i-1]:=lociq[i-1]+1;
           locik[i-1]:=1;
          END;
       END;
     end;
    newtime:=newtime+1;
    if dma then time:=time+3 else time:=time+1;
  UNTIL  time > genonum;     {genonum/3 for disease marker associate}
 likelihood:=like;
 if dg then
 if not ehplus then
  begin
   assign(tmp,'temp.dat');
   rewrite(tmp);
   write(tmp,'2 ');
   for i:=2 to totalloci do write(tmp,' ',loci[i]:2);
   writeln(tmp);
   outind(tmp,locigeno,dg11,2);
   outind(tmp,locigeno,dg12,2);
   outind(tmp,locigeno,dg22,2);
   close(tmp);
  end else
  begin
    time:=1;j:=1;
    for i:=2 to totalloci do time:=time*locigeno[i];
    for i:=1 to time do if idsave[j]=i then
    begin
      indivN[j]:=dg11[j];
      indivN[j+obscom]:=dg12[j];
      indivN[j+obscom*2]:=dg22[j];
      id[j]:=idsave[j];
      id[j+obscom]:=time+idsave[j];
      id[j+obscom*2]:=time*2+idsave[j];
      j:=j+1;
    end;
  end;
 END;   {likelihood}

 PROCEDURE gethapN;
  {compute haplotype}
 VAR i,j,time,first,last,h,subtime,subn:integer;
  linek,lineq,linep,lines:integer;
  vhaplop:linetype;

  PROCEDURE multiadd(i:integer;VAR line:integer);
   {deal with one cell which has more than one heterozygous}
  VAR j,k:integer;

  BEGIN
   IF locihn[i]=1 THEN
    BEGIN
     IF i=first THEN
      BEGIN IF i=last THEN
       BEGIN
        line:=line+1;
        {
        genotype:=true;
        linep:=linenum(locigeno,locip,1);
        genotype:=false;
        }
        linep:=saveid;
        if linep<>-1 then
        begin
          linek:=linenum(loci,locik,1);
          lineq:=linenum(loci,lociq,1);
          newhaplo[linek]:=newhaplo[linek]+vhaplop[line]*indivN[linep];
          newhaplo[lineq]:=newhaplo[lineq]+vhaplop[line]*indivN[linep];
        end;
       END
       ELSE multiadd(i+1,line);
      END
     ELSE BEGIN
      FOR j:=1 TO 2 DO
       BEGIN IF j=2 THEN swap(locik[i],lociq[i]);
        IF i=last THEN
         BEGIN
          line:=line+1;
          {
          genotype:=true;
          linep:=linenum(locigeno,locip,1);
          genotype:=false;
          }
          linep:=saveid;
          if linep<>-1 then
          begin
            linek:=linenum(loci,locik,1);
            lineq:=linenum(loci,lociq,1);
            newhaplo[linek]:=newhaplo[linek]+vhaplop[line]*indivN[linep];
            newhaplo[lineq]:=newhaplo[lineq]+vhaplop[line]*indivN[linep];
          end;
         END
        ELSE multiadd(i+1,line);
       END;
      swap(locik[i],lociq[i]);
     END;
    END
   ELSE multiadd(i+1,line);
  END;


  PROCEDURE multigen(i:integer;VAR line:integer;VAR totalg:real);
   {deal with one cell which has more than one heterozygous}
  VAR j,k:integer;

  BEGIN
   IF locihn[i]=1 THEN
    BEGIN
     IF i=first THEN
      IF i=last THEN
       BEGIN
        line:=line+1;
        linek:=linenum(loci,locik,1);
        lineq:=linenum(loci,lociq,1);
        vhaplop[line]:=2*oldhaplo[linek]*oldhaplo[lineq];
        totalg:=totalg+vhaplop[line];
       END
      ELSE multigen(i+1,line,totalg)
     ELSE BEGIN
      FOR j:=1 TO 2 DO
       BEGIN IF j=2 THEN swap(locik[i],lociq[i]);
        IF i=last THEN
         BEGIN
          line:=line+1;
          linek:=linenum(loci,locik,1);
          lineq:=linenum(loci,lociq,1);
          vhaplop[line]:=2*oldhaplo[linek]*oldhaplo[lineq];
          totalg:=totalg+vhaplop[line];
         END
        ELSE multigen(i+1,line,totalg);
       END;
      swap(locik[i],lociq[i]);
     END;
   END ELSE multigen(i+1,line,totalg);
  END;


 BEGIN  {gethapN}
  iniindex;
  newtoold;  {reset newhaplo}
  time:=1;j1:=1;j2:=1;
  REPEAT
   if ehplus then begin
      if casestudy then begin
         if time=id[j1] then begin
            saveid:=j1;j1:=j1+1;end else saveid:=-1;
      end else
      begin
        if time=idsave[j2] then begin
           saveid:=j2;j2:=j2+1;end else saveid:=-1;
      end;
      end else saveid:=time;
   for j:=1 to totalloci do
   FOR i:=1 TO totalloci DO locihn[i]:=0;
   {
   genotype:=true;
   lines:=linenum(locigeno,locip,1);
   genotype:=false;
   }
   FOR i:=1 TO totalloci DO IF (locik[i]<>lociq[i]) THEN locihn[i]:=1;
   h:=0;
   FOR i:=1 TO totalloci DO h:=h+locihn[i];
   IF h<2 THEN
    BEGIN
     {
     genotype:=true;
     linep:=linenum(locigeno,locip,1);
     genotype:=false;
     }
     linep:=saveid;
     if linep<>-1 then
     begin
       linek:=linenum(loci,locik,1);
       lineq:=linenum(loci,lociq,1);
       newhaplo[linek]:=newhaplo[linek]+indivN[linep];
       newhaplo[lineq]:=newhaplo[lineq]+indivN[linep];
     end;
    END
   ELSE IF h>1 THEN
    BEGIN
     getfirst(first,last);
     lines:=0;
     totalg:=0;
     multigen(first,lines,totalg);
     subn:=exp2(h-1);
     if totalg < verysmall then
     FOR i:=1 TO subn DO vhaplop[i]:=0
     else for i:=1 to subn do vhaplop[i]:=vhaplop[i]/totalg;
     lines:=0;
     multiadd(first,lines);
    END;
    {
    locip[totalloci]:=locip[totalloci]+1;
    FOR i:=totalloci DOWNTO 2 DO
    IF locip[i]>locigeno[i] THEN
     BEGIN
      locip[i]:=1;
      locip[i-1]:=locip[i-1]+1;
     END;
    }

  locik[totalloci]:=locik[totalloci]+1;
  if (totalloci>1) then
  begin  {for more than one marker}
   FOR i:=totalloci DOWNTO 2 DO
    IF locik[i]>lociq[i] THEN
     BEGIN
      lociq[i]:=lociq[i]+1;
      locik[i]:=1;
      IF lociq[i]>loci[i] THEN
       BEGIN
        lociq[i]:=1;
        locik[i-1]:=locik[i-1]+1;
        IF i-1=1 THEN
         IF locik[i-1]>lociq[i-1] THEN
          BEGIN lociq[i-1]:=lociq[i-1]+1;
           locik[i-1]:=1;
          END;
       END;
     END;
  end else  {for only one marker}
     IF locik[i]>lociq[i] THEN
     BEGIN
      lociq[i]:=lociq[i]+1;
      locik[i]:=1;
     end;
   time:=time+1;
  UNTIL time>genonum
 END;   {gethapN}

 FUNCTION getdiff:real;
 VAR i:integer;
  temp:real;
 BEGIN
  temp:=0;
  FOR i:=1 TO haplnum DO temp:=temp+abs(oldhaplo[i]-newhaplo[i]);
  getdiff:=temp;
 END;

procedure outdat(var outf:text);
var ss:string;
begin
  while not eof(outf) do
  begin
    readln(outf,ss);
    writeln(wf,ss);
  end;
reset(outf);
readln(outf);
end; {outdat}

procedure maxval;
begin
 writeln('The following maximum values are in effect:');
 writeln(maxalle:10,' alleles at any one locus');
 writeln(maxloci:10,' loci');
 writeln(maxgeno:10,' genotypes, ie. product over #genotypes at all loci');
 writeln(maxhap:10,' haplotypes');
 writeln;
 writeln('The program manual may be accessed at http://linkage.rockefeller.edu/ott/eh.htm');
 writeln;
end; {maxval}

FUNCTION getp(Var l:larr;k,totalloci:integer;VAR p:arr;VAR line:integer;Var ll:integer):real;
{This procedure computes alleles frequencies}
VAR i,j:integer;
    n,tempn,temp:real;

BEGIN
 tempN:=0;
 FOR i:=1 TO l[k] DO
  FOR j:=1 TO i DO
   IF k=totalloci THEN
    BEGIN
     line:=line+1;
     IF line=id[ll] then begin n:=indivn[ll];ll:=ll+1;end else n:=0;
     p[k,i]:=p[k,i]+n;
     p[k,j]:=p[k,j]+n;
     tempN:=tempN+n;
    END
   ELSE BEGIN
    temp:=getp(l,k+1,totalloci,p,line,ll);
    p[k,i]:=p[k,i]+temp;
    p[k,j]:=p[k,j]+temp;
    tempN:=tempN+temp;
   END;
 getp:=tempN;
END;

{$F+}
{$I ERRTRAP.PAS}
{$F-}

BEGIN   {eh}
 exitsave:=exitproc;
 exitproc:=@errtrap;

 cls;
 ehlogo;
 maxval;
 ch:='n';
 writeln('Do you wish to use the case-control sampling option?  [N]');
 if seekeoln then readln else readln(ch);
 casestudy:= ch in ['y','Y'];
 if ehplus then goto plus;
 if casestudy then
 begin
  fname:='control.dat';
  repeat
   done:=true;
   writeln('Enter name of control data file  [CONTROL.DAT]');
   IF seekeoln THEN readln ELSE readln(fname);
   assign(controlf,fname);
   {$I-}reset(controlf);{$I+}
   IF IOresult<>0 THEN
   begin
    done:=false;
   writeln('        File ',fname,' not found.  Try again',chr(7));
   fname:='control.dat'; end;
  until done;
   fname:='case.dat';
   repeat
     done:=true;
     writeln('Enter name of case data file   [CASE.DAT]');
     IF seekeoln THEN readln ELSE readln(fname);
     assign(casef,fname);
     {$I-}reset(casef);{$I+}
     IF IOresult<>0 THEN
     begin
       done:=false;
       writeln('        File ',fname,' not found.  Try again',chr(7));
     end;
   until done;
 end else
 begin
  fname:='eh.dat';
  repeat
   done:=true;
   writeln('Enter name of data file  [EH.DAT]');
   IF seekeoln THEN readln ELSE readln(fname);
   assign(controlf,fname);
   {$I-}reset(controlf);{$I+}
   IF IOresult<>0 THEN
   begin
    done:=false;
    writeln('        File ',fname,' not found.  Try again',chr(7));
   end;
  until done;
 end;

 writeln('Enter name of output file.  [EH.OUT]');
 fnameout:='eh.out';
 IF seekeoln THEN readln ELSE readln(fnameout);
 assign(wf,fnameout);
 {$I-}reset(wf);{$I+}
 IF NOT(IOresult<>0)
 THEN BEGIN
  writeln('File ',fnameout,' exists. Overwrite it? [y/n]',chr(7));
  IF seekeoln THEN BEGIN
   readln;
   rewrite(wf);
  END
  ELSE BEGIN
   readln(ch);
   IF ch IN ['y','Y'] THEN rewrite(wf)
   ELSE BEGIN
    programaborted;
    GOTO 999;
   END;
  END;
 END
 ELSE rewrite(wf);


{*****************start}
 cls;
 ehlogo;
 writeln('        Program is running');
 writeln;

 if casestudy then
 begin  {casestudy}

 {check control.dat}
  IF NOT inputdataok(controlf,output,totalloci,loci) THEN
  BEGIN
   programaborted;
   writeln(' Reason: Number of input genotypes should be ',genonum);
   GOTO 999;
  END ELSE BEGIN
   reset(controlf);
   readln(controlf);
  END;

 {check case.dat}
  IF NOT inputdataok(casef,output,cnloci,caseloci) THEN
  BEGIN
   programaborted;
   writeln(' Reason: Number of input genotypes should be ',genonum);
   GOTO 999;
  END ELSE BEGIN
   reset(casef);
   readln(casef);
  END;

 {check data consistant}
  if not dataconsist(loci,caseloci,totalloci,cnloci) then
  BEGIN
   programaborted;
   writeln(' Reason: Input data files are not consistent');
   GOTO 999;
  END;

  getlocih(totalloci);
  getfreq;
  i:=0;
  casecontrol(controlf,casef,loci,1,totalloci,i);
{compute initial haplotype disrecarding disease genotype}
   fname:='temp.dat';
   assign(tmp,fname);
   rewrite(tmp);
   for i:=1 to totalloci do write(tmp,' ',loci[i]:2);
   writeln(tmp);
   outind(tmp,locigeno,dg33,1);
   close(tmp);
   assign(controlf,'temp.dat');
   reset(controlf);
   IF NOT inputdataok(controlf,output,totalloci,loci) THEN
   BEGIN
    programaborted;
    writeln(' Reason: Number of input genotypes should be ',genonum);
    GOTO 999;
   END ELSE BEGIN
   reset(controlf);
   readln(controlf);
   END;
{**}
   initdg;  {write dg11,dg12 and dg22}
   getlocih(totalloci);
   initial(p);
   i:=0;
   totalind:=getN(controlf,loci,1,totalloci,p,i);
   close(controlf);

   totalall:=2*totalind;
   calcula(totalloci);  {calculate allele frequency}
   gethapP(totalloci);  {haplotype frequency based on allele frequency}

   {iniloglike:=likelihood(newhaplo,totalloci,false);}
   onemark:=totalloci=1;

   if not onemark then {allow allelic association}
   begin
     loop:=0;
     REPEAT
      loop:=loop+1;
      gethapN; {haplotype}
      FOR i:=1 TO haplnum DO newhaplo[i]:=newhaplo[i]/totalall;

      loglike:=likelihood(newhaplo,totalloci,false,false);
      diff:=getdiff;
      write('Iteration=',loop:2,'  log likelihood=',loglike:8:5,
        '   difference=',diff:8:5);
      FOR i:=1 TO 65 DO write(bs);
      done:=(loop>100) OR (diff<0.0001);  {verysmall);}
     UNTIL done;
     IF loop>100 THEN writeln('Warning: bad data set',chr(7));
   end;  {not onemark}

   for i:=1 to haplnum do
   begin
    indephaplo[i]:=newhaplo[i]*(1-pp);  {disease independent of markers}
    indephaplo[i+haplnum]:=newhaplo[i]*pp;
   end;

loop:=0;
repeat
   {/** case-control study **/}
 assign(controlf,'temp.dat');
 reset(controlf);

 IF NOT inputdataok(controlf,output,totalloci,loci) THEN
  BEGIN
   programaborted;
   writeln(' Reason: Number of input genotypes should be ',genonum);
   GOTO 999;
 END ELSE BEGIN
  reset(controlf);
  readln(controlf);
 END;

 i:=0;

 getlocih(totalloci);
 initial(p);
 totalind:=getN(controlf,loci,1,totalloci,p,i);
 close(controlf);
 totalall:=2*totalind;
 calcula(totalloci);


 if  (loop=0) then
 begin  {casestudy}
   p[1,1]:=1-pp;   {/* imporse disease gene freq. */}
   p[1,2]:=pp;

  outindP(totalloci);  {/* print allele freq. */}
  gethapP(totalloci);  {/* compute haplotype freq.s assuming no allelic asso.*/}

  iniloglike:=likelihood(newhaplo,totalloci,true,false);
  indloglike:=likelihood(indephaplo,totalloci,true,false);
  loglike:=iniloglike;
  end;

  oldlog:=loglike;
  loop1:=0;
  repeat
    gethapN;  {/* compute new haplotype freq.s */}
    FOR i:=1 TO haplnum DO newhaplo[i]:=newhaplo[i]/totalall;

    normalizeH(loci,newhaplo);
    diff:=getdiff; {put here}

    write('Iteration=',loop1:2,'  log likelihood=',loglike:8:5,'   difference=',diff:8:5);
    FOR i:=1 TO 65 DO write(bs);

    loop1:=loop1+1;
  until ((diff<0.0001) or (loop1>100));  {changed here}
  loop:=loop+1;
{  outvec(loci,newhaplo);}
  loglike:=likelihood(newhaplo,totalloci,true,true);
  diff:=loglike-oldlog;
  until ((diff<0.0001) or (loop>15));

{  outvec(loci,newhaplo);}
 end  {casestudy}
else
begin  {regular data}
 IF NOT inputdataok(controlf,output,totalloci,loci) THEN
  BEGIN
   programaborted;
   writeln(' Reason: Number of input genotypes should be ',genonum);
   GOTO 999;
 END ELSE BEGIN
  reset(controlf);
  readln(controlf);
 END;
 getlocih(totalloci);
 initial(p);
 i:=0;
 totalind:=getN(controlf,loci,1,totalloci,p,i);
 close(controlf);
 totalall:=2*totalind;
 calcula(totalloci);

  loop:=0;
  outindP(totalloci);
  gethapP(totalloci);
  iniloglike:=likelihood(newhaplo,totalloci,false,false);

  REPEAT
   loop:=loop+1;
   gethapN;
   FOR i:=1 TO haplnum DO newhaplo[i]:=newhaplo[i]/totalall;
   loglike:=likelihood(newhaplo,totalloci,false,false);
{   outvec(loci,newhaplo);}
   diff:=getdiff;
   write('Iteration=',loop:2,'  log likelihood=',loglike:8:5,'   difference=',diff:8:5);
   FOR i:=1 TO 65 DO write(bs);
  done:=(loop>100) OR (diff<0.0001);  {verysmall}
  UNTIL done;
 end;  {regular data}
 outvec(loci,newhaplo);
 close(wf);
 programsuccess;
 IF loop>100 THEN writeln('Warning: bad data set? Max. #loops exceeded',chr(7));
 goto 999;

plus:
 fname:='ehplus.dat';
 repeat
   done:=true;
   writeln('Enter name of case data file   [EHPLUS.DAT]');
   IF seekeoln THEN readln ELSE readln(fname);
   assign(controlf,fname);
   {$I-}reset(controlf);{$I+}
   IF IOresult<>0 THEN
   begin
     done:=false;
     writeln('        File ',fname,' not found.  Try again',chr(7));
   end;
 until done;

 writeln('Enter name of output file.  [EHPLUS.OUT]');
 fnameout:='ehplus.out';
 IF seekeoln THEN readln ELSE readln(fnameout);
 assign(wf,fnameout);
 {$I-}reset(wf);{$I+}
 IF NOT(IOresult<>0)
 THEN BEGIN
  writeln('File ',fnameout,' exists. Overwrite it? [y/n]',chr(7));
  IF seekeoln THEN BEGIN
   readln;
   rewrite(wf);
  END
  ELSE BEGIN
   readln(ch);
   IF ch IN ['y','Y'] THEN rewrite(wf)
   ELSE BEGIN
    programaborted;
    GOTO 999;
   END;
  END;
 END
 ELSE rewrite(wf);
 cls;
 ehlogo;
 writeln('        Program is running');
 writeln;
 numloci:=0;
 maxall:=0;
{
 for i:=1 to PRIME do
 begin
   bin[i].n:=0;
   for j:=0 to maxbinsize do begin bin[i].id[j]:=0; bin[i].obs[j]:=0; end;
 end;
}
 i:=0;
 totalind:=0;
 bb:=false;
 if casestudy then
 begin
   while not seekeof(controlf) do
   if not bb then
   begin
     while not seekeoln(controlf) do
     begin
       numloci:=numloci+1;
       read(controlf,loci[numloci]);
       if(loci[numloci]>maxall) then maxall:=loci[numloci];
     end;readln(controlf);
     bb:=true;
   end else
   begin
     i:=i+1;
     readln(controlf,idsave[i],indivN[i],caseN[i],contrlN[i]);
{
     j:=idsave[i] mod prime;
     j1:=bin[j].n;
     bin[j].obs[j1]:=i;
     bin[j].id[j1]:=idsave[i];
     bin[j].n:=bin[j].n+1;
}
     totalind:=totalind+indivN[i];id[i]:=idsave[i];
   end;
 end else
 begin
   while not seekeof(controlf) do
   if not bb then
   begin
     while not seekeoln(controlf) do
     begin
       numloci:=numloci+1;
       read(controlf,loci[numloci]);
       if(loci[numloci]>maxall) then maxall:=loci[numloci];
     end;readln(controlf);
     bb:=true;
   end else
   begin
     i:=i+1;
     readln(controlf,idsave[i],indivN[i]);
{
     j:=idsave[i] mod prime;
     j1:=bin[j].n;
     bin[j].obs[j1]:=i;
     bin[j].id[j1]:=idsave[i];
     bin[j].n:=bin[j].n+1;
}
     totalind:=totalind+indivN[i];id[i]:=idsave[i];
   end;
 end;
 close(controlf);
 for j:=1 to PRIME do
 begin
 {
   write(j:5,bin[j].n:3);
   for j1:=0 to bin[j].n do write(bin[j].id[j1]:3);writeln;
 }
 end;
 sample_size:=i;
 totalloci:=numloci;
 inivar(numloci,loci);
 obscom:=i;
 if casestudy then
 begin
   getfreq;
   for i:=1 to sample_size do
   begin
     dg11[i]:=casep[1]*caseN[i]+controlp[1]*contrlN[i];
     dg12[i]:=casep[2]*caseN[i]+controlp[2]*contrlN[i];
     dg22[i]:=casep[3]*caseN[i]+controlp[3]*contrlN[i];
   end;
   getlocih(totalloci);
   initial(p);
   i:=0;j:=1;
   totalind:=getp(loci,1,totalloci,p,i,j);
   totalall:=2*totalind;
   calcula(totalloci);
   gethapP(totalloci);
   onemark:=totalloci=1;
   if not onemark then
   begin
     loop:=0;
     REPEAT
       loop:=loop+1;
       gethapN;
       FOR i:=1 TO haplnum DO newhaplo[i]:=newhaplo[i]/totalall;
        loglike:=likelihood(newhaplo,totalloci,false,false);
        diff:=getdiff;
        write('Iteration=',loop:2,'  log likelihood=',loglike:8:5,
          '   difference=',diff:8:5);
        FOR i:=1 TO 65 DO write(bs);
        done:=(loop>100) OR (diff<0.0001);
     UNTIL done;
     IF loop>100 THEN writeln('Warning: bad data set',chr(7));
   end;
   for i:=1 to haplnum do
   begin
    indephaplo[i]:=newhaplo[i]*(1-pp);
    indephaplo[i+haplnum]:=newhaplo[i]*pp;
   end;
   j:=1;
   n:=1;
   for i:=1 to totalloci do n:=n*locigeno[i];
   for i:=totalloci+1 downto 2 do loci[i]:=loci[i-1];loci[1]:=2;
   for i:=1 to n do if(idsave[j]=i) then
   begin
     indivN[j]:=dg11[j];
     indivN[obscom+j]:=dg12[j];
     indivN[obscom*2+j]:=dg22[j];
     id[j]:=idsave[j];
     id[obscom+j]:=n+idsave[j];
     id[obscom*2+j]:=2*n+idsave[j];
{
     j1:=id[j] mod prime;
     j2:=bin[j1].n;
     bin[j1].obs[j2]:=j;
     bin[j1].id[j2]:=id[j];
     bin[j1].n:=bin[j1].n+1;
     j1:=id[obscom+j] mod prime;
     j2:=bin[j1].n;
     bin[j1].obs[j2]:=obscom+j;
     bin[j1].id[j2]:=id[obscom+j];
     bin[j1].n:=bin[j1].n+1;
     j1:=id[obscom*2+j] mod prime;
     j2:=bin[j1].n;
     bin[j1].obs[j2]:=obscom*2+j;
     bin[j1].id[j2]:=id[obscom*2+j];
     bin[j1].n:=bin[j1].n+1;
}
     j:=j+1;
   end;
   totalloci:=totalloci+1;
   numloci:=numloci+1;
   haplnum:=2*haplnum;
   inivar(numloci,loci);
   loop:=0;
   repeat
     getlocih(totalloci);
     initial(p);
     i:=0;j:=1;
     totalind:=getp(loci,1,totalloci,p,i,j);
     totalall:=2*totalind;
     calcula(totalloci);
     if(loop=0) then
     begin
       p[1,1]:=1-pp;
       p[1,2]:=pp;
       outindP(totalloci);
       gethapP(totalloci);
       iniloglike:=likelihood(newhaplo,totalloci,true,false);
       indloglike:=likelihood(indephaplo,totalloci,true,false);
       loglike:=iniloglike;
     end;
     oldlog:=loglike;
     loop1:=0;
     repeat
       gethapN;
       FOR i:=1 TO haplnum DO newhaplo[i]:=newhaplo[i]/totalall;
       normalizeH(loci,newhaplo);
       diff:=getdiff;
       write('Iteration=',loop1:2,'  log likelihood=',loglike:8:5,'   difference=',diff:8:5);
       FOR i:=1 TO 65 DO write(bs);
       loop1:=loop1+1;
     until ((diff<0.0001) or (loop1>100));
     loop:=loop+1;
     loglike:=likelihood(newhaplo,totalloci,true,true);
     diff:=loglike-oldlog;
   until ((diff<0.0001) or (loop>15));
 end else
 begin
   getlocih(totalloci);
   initial(p);
   i:=0;j:=1;
   totalind:=getp(loci,1,totalloci,p,i,j);
   totalall:=2*totalind;
   calcula(totalloci);
   loop:=0;
   outindP(totalloci);
   gethapP(totalloci);
   iniloglike:=likelihood(newhaplo,totalloci,false,false);
   REPEAT
     loop:=loop+1;
     gethapN;
     FOR i:=1 TO haplnum DO newhaplo[i]:=newhaplo[i]/totalall;
     loglike:=likelihood(newhaplo,totalloci,false,false);
     diff:=getdiff;
     write('Iteration=',loop:2,'  log likelihood=',loglike:8:5,'   difference=',diff:8:5);
     FOR i:=1 TO 65 DO write(bs);
     done:=(loop>100) OR (diff<0.0001);
   UNTIL done;
 end;
 outvec(loci,newhaplo);
 close(wf);
 programsuccess;
 IF loop>100 THEN writeln('Warning: bad data set? Max. #loops exceeded',chr(7));
999:
END. {eh}