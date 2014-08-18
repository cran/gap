/*9-12-9 MRC-Epid JHZ*/

proc proto package = work.funcs label = "HWE for SNPs";
     #define NULL 0;
     double SNPHWE(int obs_hets, int obs_hom1, int obs_hom2);

     externc SNPHWE;
     double SNPHWE(int obs_hets, int obs_hom1, int obs_hom2)
     {
     int obs_homc, obs_homr, rare_copies, genotypes;
     int i, mid;
     int curr_hets, curr_homr, curr_homc;
     double sum, p_hwe, *het_probs;
     if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0) p_hwe = -1;
     obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
     obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;
     rare_copies = 2 * obs_homr + obs_hets;
     genotypes   = obs_hets + obs_homc + obs_homr;
     het_probs = (double *)malloc((sizeof(double)*(rare_copies + 1)));
     if (het_probs == NULL) p_hwe = -1;
     for (i = 0; i <= rare_copies; i++) het_probs[i] = 0.0;
  /* start at midpoint */
     mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);
  /* check to ensure that midpoint and rare alleles have same parity */
     if ((rare_copies & 1) ^ (mid & 1)) mid++;
     curr_hets = mid;
     curr_homr = (rare_copies - mid) / 2;
     curr_homc = genotypes - curr_hets - curr_homr;
     het_probs[mid] = 1.0;
     sum = het_probs[mid];
     for (curr_hets = mid; curr_hets > 1; curr_hets -= 2)
     {
         het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0)
                               / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
         sum += het_probs[curr_hets - 2];
      /* 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote */
         curr_homr++;
         curr_homc++;
     }
     curr_hets = mid;
     curr_homr = (rare_copies - mid) / 2;
     curr_homc = genotypes - curr_hets - curr_homr;
     for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2)
     {
         het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc
                            /((curr_hets + 2.0) * (curr_hets + 1.0));
         sum += het_probs[curr_hets + 2];
      /* add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote */
         curr_homr--;
         curr_homc--;
     }
     for (i = 0; i <= rare_copies; i++)  het_probs[i] /= sum;
  /* alternate p-value calculation for p_hi/p_lo
     double p_hi = het_probs[obs_hets];
     for (i = obs_hets + 1; i <= rare_copies; i++)
     p_hi += het_probs[i];
     double p_lo = het_probs[obs_hets];
     for (i = obs_hets - 1; i >= 0; i--) p_lo += het_probs[i];
     double p_hi_lo = p_hi < p_lo ? 2.0 * p_hi : 2.0 * p_lo;
  */
     p_hwe = 0.0;
     for (i = 0; i <= rare_copies; i++)
     {
         if (het_probs[i] > het_probs[obs_hets]) continue;
         p_hwe += het_probs[i];
     }
     p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;
     free(het_probs);
     return p_hwe;
     }
     externcend;
run;

proc fcmp inlib=work outlib=work.funcs.trial;
     function HWE(b,a,c);
              pHWE=SNPHWE(b,a,c);
              return(pHWE);
     endsub;
     result = SNPHWE(2,1,3);
     put result=;
     result = SNPHWE(20,10,30);
     put result=;
     result = SNPHWE(200,100,300);
     put result=;
run;

options cmplib=(work work.funcs);
data abc;
     input a b c;
     pHWE=HWE(b,a,c);
     datalines;
     1 2 3
     10 20 30
     100 200 300
run;

options nocenter;
proc print;
     format pHWE 20.15;
run;
