\name{pedtodot}
\alias{pedtodot}
\title{Converting pedigree(s) to dot file(s)}
\description{
This function converts GAS or LINKAGE formatted pedigree(s) into .dot file
for each pedigree to be used by dot in graphviz, which is a flexible package
for graphics freely available from \url{http://www.graphviz.org}

Note that a single PostScript file is obtainable by specifyiing *.dot to dot.
However, to preserve the original order of pedigree(s) in the data, you can
examine the examples at the end of this document.

Under Cygwin/Linux/Unix, the PostScript file can be converted to Portable 
Document Format (PDF) default to Acrobat.

\% ps2pdf <ps file>

Use ps2pdf12, ps2pdf13, or ps2pdf14 for appropriate versions of Acrobat
according to information given on the headline of <ps file>.
}

\usage{
pedtodot(pedfile,makeped=F,sink=T,page="",url="",height=0.5,width=0.75,rotate=90,dir="none")
}

\arguments{
\item{pedfile}{a pedigree file in GAS or LINKAGE format, note if 
individual's ID is character then it is necessary to specify as.is=T
in the read.table command}
\item{makeped}{a logical variable indicating if the pedigree file is post-makeped}
\item{sink}{a logical variable indicating if .dot file(s) are created}
\item{page}{a string indicating the page size, e.g, A4, A5, B5, Legal, Letter, 
Executive, "x,y", where x, y is the customized page size}
\item{url}{Unified Resource Locator (URL) associated with the diagram(s)}
\item{height}{the height of node(s)}
\item{width}{the width of node(s)}
\item{rotate}{if set to 90, the diagram is in landscape}
\item{dir}{direction of edges, i.e., "none", "forward","back","both". This will be useful
if the diagram is viewed by lneato}
}

\details{
We can extract the code below (or within pedtodot.Rd) to pedtodot and then
use command: 

sh pedtodot <pedigree file>

\preformatted{
# Read a GAS or LINKAGE format pedigree, return a digraph in the dot language 
# call dot to make pedigree drawing
# 
AWK=/usr/local/bin/gawk
DOTEXE=/people/jzhao/gv1.7c/bin/dot
# cygwin
# AWK=/bin/gawk
# DOTEXE=c:/local/graphviz/bin/dot

for fil in $*
do
  for ped in `$AWK '!/^[!#]/ {print $1}' $fil | sort -u`
  do
     echo "Pedigree $ped"
     $AWK -v ped=$ped '
     BEGIN { shape["m"]="box,regular=1"
             shape["1"]="box,regular=1"
             shape["f"]="circle"
             shape["2"]="circle"
             shade["y"]="blue"
             shade["2"]="blue"
             shade["n"]="grey"
             shade["1"]="grey"
             shade["x"]="green"
             shade["0"]="green"
     }
     !/^[!#]/ && $1==ped {
             sex[$2]=$5
             aff[$2]="x" ; if ($6 ~ /[012nyx]/) aff[$2]=$6
             if($3!="x" && $3!="0") {
               marriage[$3,$4]++
               child[$3,$4,marriage[$3,$4]]=$2
             }
     }
     END   { print "digraph Ped_" ped " {"
             print "# page =\"8.2677165,11.692913\" ;"
             print "ratio =\"auto\" ;"
             print "mincross = 2.0 ;"
             print "label=\"Pedigree " ped "\" ;"
             print "rotate=90 ;"
             for(s in sex) {
               print "\"" s "\" [shape=" shape[sex[s]] ","  \
                     " style=filled,color=" shade[aff[s]] "] ;"
             }
             for(m in marriage) {
               n=split(m,par,"\034")
               mating="\"" par[1] "x" par[2] "\""
               print mating "[shape=diamond,style=filled," \
                     "label=\"\",height=.1,width=.1] ;"
               print "\"" par[1] "\" -> " mating " [dir=none, weight=1] ;"
               print "\"" par[2] "\" -> " mating " [dir=none, weight=1] ;"
               for(k=1;k<=marriage[par[1],par[2]];k++) {
                 print  mating " -> \"" child[par[1],par[2],k] "\"" \
                        " [dir=none, weight=2] ;"
               }
             }
             print "}"
     }' $fil > $ped.dot 
     $DOTEXE -Tps $ped.dot -o $ped.ps
  done
done
$DOTEXE -Tps *.dot -o $*.ps
}
}

\value{
For each pedigree, the function generates a .dot file to be used by dot. The
collection of all pedigrees (*.dot) can also be put together.
}

\seealso{package sem in CRAN and Rgraphviz in BioConductor \url{http://www.bioconductor.org}}

\examples{\dontrun{
#
# Three examples of pedigree-drawing
# assuming pre-MakePed LINKAGE file in which IDs are characters
pre<-read.table("pheno.pre",as.is=T)[,1:6]
pedtodot(pre)
dir()      
# for post-MakePed LINKAGE file in which IDs are integers
ped <-read.table("pheno.ped")[,1:10]
pedtodot(ped,makeped=T)
dir()
# for a single file with a list of pedigrees ordered data
sink("gaw14.dot")
pedtodot(ped,sink=F)
sink()
file.show("gaw14.dot")
# more details
pedtodot(ped,sink=F,page="B5",url="http://www.ucl.ac.uk/~rmjdjhz")
}
}
\author{David Duffy, Jing Hua Zhao}
\note{This is based on the gawk script program pedtodot by David Duffy with minor changes}
\keyword{dplot}
