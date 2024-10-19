* Stata codes for figures and endogenous kink model publication bias test

* Import data 
clear
import excel "C:\Users\USER PC\Documents\Data.xlsx", sheet("All") firstrow


* Codes for funnel plots and box plots

twoway (scatter inv_se beta, mcolor(navy) msymbol(smcircle_hollow)) if hike == 1 & interbank ==1, ytitle(`"Precision: 1/standard error"') xtitle(`"Beta: percentage increase in lending rate due to a 1% policy rate hike"') xline(0.6841, lpattern(solid) lcolor(red) lwidth(vthin)) xline(0.5530, lpattern(dash) lcolor(red) lwidth(vthin)) title(`"Upward pass-through: interbank rates"') xsize(6.5) ysize(5.0)


twoway (scatter inv_se beta, mcolor(navy) msymbol(smcircle_hollow)) if cut == 1 & interbank ==1, ytitle(`"Precision: 1/standard error"') xtitle(`"Beta: percentage increase in lending rate due to a 1% policy rate cut"') xline(0.1206, lpattern(solid) lcolor(red) lwidth(vthin)) xline(0.1152, lpattern(dash) lcolor(red) lwidth(vthin)) title(`"Downward pass-through: interbank rates"') xsize(6.5) ysize(5.0)



twoway (scatter inv_se beta, mcolor(navy) msymbol(smcircle_hollow)) if hike == 1 & discount ==1, ytitle(`"Precision: 1/standard error"') xtitle(`"Beta: percentage increase in lending rate due to a 1% policy rate hike"') xline(0.0645, lpattern(solid) lcolor(red) lwidth(vthin)) xline(0.0130, lpattern(dash) lcolor(red) lwidth(vthin)) title(`"Upward pass-through: discount rates"') xsize(6.5) ysize(5.0)



twoway (scatter inv_se beta, mcolor(navy) msymbol(smcircle_hollow)) if cut == 1 & discount ==1, ytitle(`"Precision: 1/standard error"') xtitle(`"Beta: percentage increase in lending rate due to a 1% policy rate cut"') xline(-0.0757, lpattern(solid) lcolor(red) lwidth(vthin)) xline(-0.0480, lpattern(dash) lcolor(red) lwidth(vthin)) title(`"Downward pass-through: discount rates"') xsize(6.5) ysize(5.0)



graph hbox beta if hike == 1, over(Study) xsize(9) ysize(11) scale(0.5)


graph hbox beta if cut == 1, over(Study) xsize(9) ysize(11) scale(0.5)






* Codes for endogenous kink model (Bom and Rachinger, 2019):

* Upward pass-through:
quietly {
	
* Refresh and import data again
clear
import excel "C:\Users\USER PC\Documents\Data.xlsx", sheet("All") firstrow


set more off
xtset study_id


* Select upward pass-through estimates

drop if cut ==1



rename beta bs
rename se sebs



gen ones=1

sum
local M=r(N)
sum sebs
local sebs_min=r(min)
local sebs_max=r(max)

gen sebs2=sebs^2
gen wis=ones/sebs2
gen bs_sebs=bs/sebs
gen ones_sebs=ones/sebs

gen bswis=bs*wis
sum wis
local wis_sum=r(sum)


* FAT-PET
regress bs_sebs ones_sebs ones,noc
local pet=_b[ones_sebs]
local t1_linreg = (_b[ones_sebs]/_se[ones_sebs])
local b_lin=_b[ones_sebs]
local Q1_lin = e(rss)
di `t1_linreg'
local abs_t1_linreg = abs(`t1_linreg')
di `abs_t1_linreg'


* PEESE
regress bs_sebs ones_sebs sebs,noc
local peese=_b[ones_sebs]
local b_sq=_b[ones_sebs]
local Q1_sq = e(rss)
di `Q1_sq'


* FAT-PET-PEESE

if `abs_t1_linreg' > invt(`M-2', 0.975) {
local combreg=`b_sq'
local Q1=`Q1_sq'
}
else {
local combreg=`b_lin'
local Q1=`Q1_lin'
}

* estimation of random effects variance component
local sigh2hat=max(0,`M'*((`Q1'/(`M'-e(df_m)-1))-1)/`wis_sum') 
local sighhat=sqrt(`sigh2hat') 


* Cutoff value for EK
if `combreg'>1.96*`sighhat' {
local a1=(`combreg'-1.96*`sighhat')*(`combreg'+1.96*`sighhat')/(2*1.96*`combreg')
}
else {
local a1=0
}



gen bs_instru = bs/se_instru
gen ones_instru = ones/se_instru
gen sebs_instru = sebs/se_instru

rename bs bs_original
rename bs_instru bs
rename ones_instru constant
rename sebs_instru pub_bias


noisily: display "EK regression: "

if `a1'>`sebs_min' & `a1'<`sebs_max' {
gen sebs_a1=sebs-`a1' if sebs>`a1'
replace sebs_a1=0 if sebs<=`a1'
gen pubbias=sebs_a1/se_instru
noisily regress bs constant pubbias, noc vce(cluster study_id)
local b0_ek=_b[constant]
local b1_ek=_b[pubbias]
local sd0_ek=_se[constant]
local sd1_ek=_se[pubbias]
}
else if `a1'<`sebs_min' {
noisily regress bs constant pub_bias, noc vce(cluster study_id)
local b0_ek=_b[constant]
local b1_ek=_b[pub_bias]
local sd0_ek=_se[constant]
local sd1_ek=_se[pub_bias]
}
else if `a1'>`sebs_max' {
noisily regress bs constant, noc vce(cluster study_id)
local b0_ek=_b[constant]
local sd0_ek=_se[constant]	
}
noisily: display "EK's mean effect estimate (alpha1) and standard error:"
noisily: di `b0_ek' 
noisily: di `sd0_ek' 
noisily: display "EK's publication bias estimate (delta) and standard error:"
noisily: di `b1_ek' 
noisily: di `sd1_ek' 


clear
}


* Downward pass-through:
quietly {
	
* Refresh and import data again
clear
import excel "C:\Users\USER PC\Documents\Data.xlsx", sheet("All") firstrow


set more off
xtset study_id


* Select downward pass-through estimates

drop if hike ==1



rename beta bs
rename se sebs



gen ones=1

sum
local M=r(N)
sum sebs
local sebs_min=r(min)
local sebs_max=r(max)

gen sebs2=sebs^2
gen wis=ones/sebs2
gen bs_sebs=bs/sebs
gen ones_sebs=ones/sebs

gen bswis=bs*wis
sum wis
local wis_sum=r(sum)


* FAT-PET
regress bs_sebs ones_sebs ones,noc
local pet=_b[ones_sebs]
local t1_linreg = (_b[ones_sebs]/_se[ones_sebs])
local b_lin=_b[ones_sebs]
local Q1_lin = e(rss)
di `t1_linreg'
local abs_t1_linreg = abs(`t1_linreg')
di `abs_t1_linreg'


* PEESE
regress bs_sebs ones_sebs sebs,noc
local peese=_b[ones_sebs]
local b_sq=_b[ones_sebs]
local Q1_sq = e(rss)
di `Q1_sq'


* FAT-PET-PEESE

if `abs_t1_linreg' > invt(`M-2', 0.975) {
local combreg=`b_sq'
local Q1=`Q1_sq'
}
else {
local combreg=`b_lin'
local Q1=`Q1_lin'
}

* estimation of random effects variance component
local sigh2hat=max(0,`M'*((`Q1'/(`M'-e(df_m)-1))-1)/`wis_sum') 
local sighhat=sqrt(`sigh2hat') 


* Cutoff value for EK
if `combreg'>1.96*`sighhat' {
local a1=(`combreg'-1.96*`sighhat')*(`combreg'+1.96*`sighhat')/(2*1.96*`combreg')
}
else {
local a1=0
}



gen bs_instru = bs/se_instru
gen ones_instru = ones/se_instru
gen sebs_instru = sebs/se_instru

rename bs bs_original
rename bs_instru bs
rename ones_instru constant
rename sebs_instru pub_bias


noisily: display "EK regression: "

if `a1'>`sebs_min' & `a1'<`sebs_max' {
gen sebs_a1=sebs-`a1' if sebs>`a1'
replace sebs_a1=0 if sebs<=`a1'
gen pubbias=sebs_a1/se_instru
noisily regress bs constant pubbias, noc vce(cluster study_id)
local b0_ek=_b[constant]
local b1_ek=_b[pubbias]
local sd0_ek=_se[constant]
local sd1_ek=_se[pubbias]
}
else if `a1'<`sebs_min' {
noisily regress bs constant pub_bias, noc vce(cluster study_id)
local b0_ek=_b[constant]
local b1_ek=_b[pub_bias]
local sd0_ek=_se[constant]
local sd1_ek=_se[pub_bias]
}
else if `a1'>`sebs_max' {
noisily regress bs constant, noc vce(cluster study_id)
local b0_ek=_b[constant]
local sd0_ek=_se[constant]	
}
noisily: display "EK's mean effect estimate (alpha1) and standard error:"
noisily: di `b0_ek' 
noisily: di `sd0_ek' 
noisily: display "EK's publication bias estimate (delta) and standard error:"
noisily: di `b1_ek' 
noisily: di `sd1_ek' 


clear
}



* Upward pass-through: interbank
quietly {
	
* Refresh and import data again
clear
import excel "C:\Users\USER PC\Documents\Data.xlsx", sheet("All") firstrow


set more off
xtset study_id


* Select upward pass-through estimates

drop if cut ==1 | discount == 1



rename beta bs
rename se sebs



gen ones=1

sum
local M=r(N)
sum sebs
local sebs_min=r(min)
local sebs_max=r(max)

gen sebs2=sebs^2
gen wis=ones/sebs2
gen bs_sebs=bs/sebs
gen ones_sebs=ones/sebs

gen bswis=bs*wis
sum wis
local wis_sum=r(sum)


* FAT-PET
regress bs_sebs ones_sebs ones,noc
local pet=_b[ones_sebs]
local t1_linreg = (_b[ones_sebs]/_se[ones_sebs])
local b_lin=_b[ones_sebs]
local Q1_lin = e(rss)
di `t1_linreg'
local abs_t1_linreg = abs(`t1_linreg')
di `abs_t1_linreg'


* PEESE
regress bs_sebs ones_sebs sebs,noc
local peese=_b[ones_sebs]
local b_sq=_b[ones_sebs]
local Q1_sq = e(rss)
di `Q1_sq'


* FAT-PET-PEESE

if `abs_t1_linreg' > invt(`M-2', 0.975) {
local combreg=`b_sq'
local Q1=`Q1_sq'
}
else {
local combreg=`b_lin'
local Q1=`Q1_lin'
}

* estimation of random effects variance component
local sigh2hat=max(0,`M'*((`Q1'/(`M'-e(df_m)-1))-1)/`wis_sum') 
local sighhat=sqrt(`sigh2hat') 


* Cutoff value for EK
if `combreg'>1.96*`sighhat' {
local a1=(`combreg'-1.96*`sighhat')*(`combreg'+1.96*`sighhat')/(2*1.96*`combreg')
}
else {
local a1=0
}



gen bs_instru = bs/se_instru
gen ones_instru = ones/se_instru
gen sebs_instru = sebs/se_instru

rename bs bs_original
rename bs_instru bs
rename ones_instru constant
rename sebs_instru pub_bias


noisily: display "EK regression: "

if `a1'>`sebs_min' & `a1'<`sebs_max' {
gen sebs_a1=sebs-`a1' if sebs>`a1'
replace sebs_a1=0 if sebs<=`a1'
gen pubbias=sebs_a1/se_instru
noisily regress bs constant pubbias, noc vce(cluster study_id)
local b0_ek=_b[constant]
local b1_ek=_b[pubbias]
local sd0_ek=_se[constant]
local sd1_ek=_se[pubbias]
}
else if `a1'<`sebs_min' {
noisily regress bs constant pub_bias, noc vce(cluster study_id)
local b0_ek=_b[constant]
local b1_ek=_b[pub_bias]
local sd0_ek=_se[constant]
local sd1_ek=_se[pub_bias]
}
else if `a1'>`sebs_max' {
noisily regress bs constant, noc vce(cluster study_id)
local b0_ek=_b[constant]
local sd0_ek=_se[constant]	
}
noisily: display "EK's mean effect estimate (alpha1) and standard error:"
noisily: di `b0_ek' 
noisily: di `sd0_ek' 
noisily: display "EK's publication bias estimate (delta) and standard error:"
noisily: di `b1_ek' 
noisily: di `sd1_ek' 


clear
}



* Downward pass-through: interbank
quietly {
	
* Refresh and import data again
clear
import excel "C:\Users\USER PC\Documents\Data.xlsx", sheet("All") firstrow


set more off
xtset study_id


* Select downward pass-through estimates

drop if hike ==1 | discount == 1



rename beta bs
rename se sebs



gen ones=1

sum
local M=r(N)
sum sebs
local sebs_min=r(min)
local sebs_max=r(max)

gen sebs2=sebs^2
gen wis=ones/sebs2
gen bs_sebs=bs/sebs
gen ones_sebs=ones/sebs

gen bswis=bs*wis
sum wis
local wis_sum=r(sum)


* FAT-PET
regress bs_sebs ones_sebs ones,noc
local pet=_b[ones_sebs]
local t1_linreg = (_b[ones_sebs]/_se[ones_sebs])
local b_lin=_b[ones_sebs]
local Q1_lin = e(rss)
di `t1_linreg'
local abs_t1_linreg = abs(`t1_linreg')
di `abs_t1_linreg'


* PEESE
regress bs_sebs ones_sebs sebs,noc
local peese=_b[ones_sebs]
local b_sq=_b[ones_sebs]
local Q1_sq = e(rss)
di `Q1_sq'


* FAT-PET-PEESE

if `abs_t1_linreg' > invt(`M-2', 0.975) {
local combreg=`b_sq'
local Q1=`Q1_sq'
}
else {
local combreg=`b_lin'
local Q1=`Q1_lin'
}

* estimation of random effects variance component
local sigh2hat=max(0,`M'*((`Q1'/(`M'-e(df_m)-1))-1)/`wis_sum') 
local sighhat=sqrt(`sigh2hat') 


* Cutoff value for EK
if `combreg'>1.96*`sighhat' {
local a1=(`combreg'-1.96*`sighhat')*(`combreg'+1.96*`sighhat')/(2*1.96*`combreg')
}
else {
local a1=0
}



gen bs_instru = bs/se_instru
gen ones_instru = ones/se_instru
gen sebs_instru = sebs/se_instru

rename bs bs_original
rename bs_instru bs
rename ones_instru constant
rename sebs_instru pub_bias


noisily: display "EK regression: "

if `a1'>`sebs_min' & `a1'<`sebs_max' {
gen sebs_a1=sebs-`a1' if sebs>`a1'
replace sebs_a1=0 if sebs<=`a1'
gen pubbias=sebs_a1/se_instru
noisily regress bs constant pubbias, noc vce(cluster study_id)
local b0_ek=_b[constant]
local b1_ek=_b[pubbias]
local sd0_ek=_se[constant]
local sd1_ek=_se[pubbias]
}
else if `a1'<`sebs_min' {
noisily regress bs constant pub_bias, noc vce(cluster study_id)
local b0_ek=_b[constant]
local b1_ek=_b[pub_bias]
local sd0_ek=_se[constant]
local sd1_ek=_se[pub_bias]
}
else if `a1'>`sebs_max' {
noisily regress bs constant, noc vce(cluster study_id)
local b0_ek=_b[constant]
local sd0_ek=_se[constant]	
}
noisily: display "EK's mean effect estimate (alpha1) and standard error:"
noisily: di `b0_ek' 
noisily: di `sd0_ek' 
noisily: display "EK's publication bias estimate (delta) and standard error:"
noisily: di `b1_ek' 
noisily: di `sd1_ek' 


clear
}



* Upward pass-through: discount
quietly {
	
* Refresh and import data again
clear
import excel "C:\Users\USER PC\Documents\Data.xlsx", sheet("All") firstrow


set more off
xtset study_id


* Select upward pass-through estimates

drop if cut ==1 | interbank == 1



rename beta bs
rename se sebs



gen ones=1

sum
local M=r(N)
sum sebs
local sebs_min=r(min)
local sebs_max=r(max)

gen sebs2=sebs^2
gen wis=ones/sebs2
gen bs_sebs=bs/sebs
gen ones_sebs=ones/sebs

gen bswis=bs*wis
sum wis
local wis_sum=r(sum)


* FAT-PET
regress bs_sebs ones_sebs ones,noc
local pet=_b[ones_sebs]
local t1_linreg = (_b[ones_sebs]/_se[ones_sebs])
local b_lin=_b[ones_sebs]
local Q1_lin = e(rss)
di `t1_linreg'
local abs_t1_linreg = abs(`t1_linreg')
di `abs_t1_linreg'


* PEESE
regress bs_sebs ones_sebs sebs,noc
local peese=_b[ones_sebs]
local b_sq=_b[ones_sebs]
local Q1_sq = e(rss)
di `Q1_sq'


* FAT-PET-PEESE

if `abs_t1_linreg' > invt(`M-2', 0.975) {
local combreg=`b_sq'
local Q1=`Q1_sq'
}
else {
local combreg=`b_lin'
local Q1=`Q1_lin'
}

* estimation of random effects variance component
local sigh2hat=max(0,`M'*((`Q1'/(`M'-e(df_m)-1))-1)/`wis_sum') 
local sighhat=sqrt(`sigh2hat') 


* Cutoff value for EK
if `combreg'>1.96*`sighhat' {
local a1=(`combreg'-1.96*`sighhat')*(`combreg'+1.96*`sighhat')/(2*1.96*`combreg')
}
else {
local a1=0
}



gen bs_instru = bs/se_instru
gen ones_instru = ones/se_instru
gen sebs_instru = sebs/se_instru

rename bs bs_original
rename bs_instru bs
rename ones_instru constant
rename sebs_instru pub_bias


noisily: display "EK regression: "

if `a1'>`sebs_min' & `a1'<`sebs_max' {
gen sebs_a1=sebs-`a1' if sebs>`a1'
replace sebs_a1=0 if sebs<=`a1'
gen pubbias=sebs_a1/se_instru
noisily regress bs constant pubbias, noc vce(cluster study_id)
local b0_ek=_b[constant]
local b1_ek=_b[pubbias]
local sd0_ek=_se[constant]
local sd1_ek=_se[pubbias]
}
else if `a1'<`sebs_min' {
noisily regress bs constant pub_bias, noc vce(cluster study_id)
local b0_ek=_b[constant]
local b1_ek=_b[pub_bias]
local sd0_ek=_se[constant]
local sd1_ek=_se[pub_bias]
}
else if `a1'>`sebs_max' {
noisily regress bs constant, noc vce(cluster study_id)
local b0_ek=_b[constant]
local sd0_ek=_se[constant]	
}
noisily: display "EK's mean effect estimate (alpha1) and standard error:"
noisily: di `b0_ek' 
noisily: di `sd0_ek' 
noisily: display "EK's publication bias estimate (delta) and standard error:"
noisily: di `b1_ek' 
noisily: di `sd1_ek' 


clear
}



* Downward pass-through: discount
quietly {
	
* Refresh and import data again
clear
import excel "C:\Users\USER PC\Documents\Data.xlsx", sheet("All") firstrow


set more off
xtset study_id


* Select downward pass-through estimates

drop if hike ==1 | interbank == 1



rename beta bs
rename se sebs



gen ones=1

sum
local M=r(N)
sum sebs
local sebs_min=r(min)
local sebs_max=r(max)

gen sebs2=sebs^2
gen wis=ones/sebs2
gen bs_sebs=bs/sebs
gen ones_sebs=ones/sebs

gen bswis=bs*wis
sum wis
local wis_sum=r(sum)


* FAT-PET
regress bs_sebs ones_sebs ones,noc
local pet=_b[ones_sebs]
local t1_linreg = (_b[ones_sebs]/_se[ones_sebs])
local b_lin=_b[ones_sebs]
local Q1_lin = e(rss)
di `t1_linreg'
local abs_t1_linreg = abs(`t1_linreg')
di `abs_t1_linreg'


* PEESE
regress bs_sebs ones_sebs sebs,noc
local peese=_b[ones_sebs]
local b_sq=_b[ones_sebs]
local Q1_sq = e(rss)
di `Q1_sq'


* FAT-PET-PEESE

if `abs_t1_linreg' > invt(`M-2', 0.975) {
local combreg=`b_sq'
local Q1=`Q1_sq'
}
else {
local combreg=`b_lin'
local Q1=`Q1_lin'
}

* estimation of random effects variance component
local sigh2hat=max(0,`M'*((`Q1'/(`M'-e(df_m)-1))-1)/`wis_sum') 
local sighhat=sqrt(`sigh2hat') 


* Cutoff value for EK
if `combreg'>1.96*`sighhat' {
local a1=(`combreg'-1.96*`sighhat')*(`combreg'+1.96*`sighhat')/(2*1.96*`combreg')
}
else {
local a1=0
}



gen bs_instru = bs/se_instru
gen ones_instru = ones/se_instru
gen sebs_instru = sebs/se_instru

rename bs bs_original
rename bs_instru bs
rename ones_instru constant
rename sebs_instru pub_bias


noisily: display "EK regression: "

if `a1'>`sebs_min' & `a1'<`sebs_max' {
gen sebs_a1=sebs-`a1' if sebs>`a1'
replace sebs_a1=0 if sebs<=`a1'
gen pubbias=sebs_a1/se_instru
noisily regress bs constant pubbias, noc vce(cluster study_id)
local b0_ek=_b[constant]
local b1_ek=_b[pubbias]
local sd0_ek=_se[constant]
local sd1_ek=_se[pubbias]
}
else if `a1'<`sebs_min' {
noisily regress bs constant pub_bias, noc vce(cluster study_id)
local b0_ek=_b[constant]
local b1_ek=_b[pub_bias]
local sd0_ek=_se[constant]
local sd1_ek=_se[pub_bias]
}
else if `a1'>`sebs_max' {
noisily regress bs constant, noc vce(cluster study_id)
local b0_ek=_b[constant]
local sd0_ek=_se[constant]	
}
noisily: display "EK's mean effect estimate (alpha1) and standard error:"
noisily: di `b0_ek' 
noisily: di `sd0_ek' 
noisily: display "EK's publication bias estimate (delta) and standard error:"
noisily: di `b1_ek' 
noisily: di `sd1_ek' 


clear
}

