cd "P:\AddHealth\Contract\27110801-Conley\Work\Luyin Zhang\admixture"
use interm_pheno.dta, clear
logit dna i.racec i.sex byear i.edu i.strataid [pweight = wt4], or
scalar loglik = e(ll)
scalar pr2 = e(r2_p)
predict prob, pr
gen ipwgt = wt4 / prob
margins, dydx(*) post
estadd scalar loglik = loglik
estadd scalar pr2 = pr2
est store m
replace ipwgt = . if dna==0
esttab m using "tables/appendix_ipw.csv", replace nogaps compress label ///
  star(* 0.05 ** 0.01 *** 0.001) b(%20.3f) se(%20.3f) obslast eqlabel(none) ///
  scalars("pr2 Pseudo R2") coeflabels(byear "Birth year (wave 4)") nonumbers mtitles("AME")
save all_pheno.dta, replace

