************************************************************
* Project 1: Verification of the Analyst Report
* Cai and Szeidl (2024), "Indirect Effects of Access to Finance"
*
* Starts from the analyst's LoanAnalysis.do. For each section:
*   (a) reproduces the analyst's result, then
*   (b) runs the corrected analysis.
* Run from the folder that contains loanmain.dta:
*   cd "<folder with loanmain.dta>"
************************************************************

clear all
set more off
capture log close
log using "LoanVerification.log", replace text

capture which esttab
if _rc ssc install estout

use "loanmain.dta", clear

* Variable definitions used below
*   type             = firm was randomly assigned the treatment (0/1)
*   survey_town_type = MARKET treatment intensity (0=pure control, 1=50%, 2=80%)
*   T50 C50 T80 C80  = treated / untreated firm in a 50% / 80% market
*                      (omitted group = pure control: untreated, in a 0% market)
*   newloan          = firm borrowed through the new program (a CHOICE, not randomized)
*   round            = 1 baseline (2013), 2 midline (2015), 3 endline (2016)
*   treatratio_comp  = share of the firm's competitors that were treated

************************************************************
* SECTION 1. DATA AND ESTIMATION SETUP
************************************************************

* Structure of the data
describe, short
isid firmid round                         // firm x round panel
tab round                                 // 3,173 firms observed in 3 rounds
bysort firmid: assert type==type[1]       // treatment is fixed within firm
bysort firmid: assert newloan==newloan[1] // newloan = ever borrowed by endline

* Market-level randomization: observations are NOT independent
egen tag_market = tag(survey_town)
count if tag_market==1                    // 78 markets
tab survey_town_type type if round==1     // pure-control markets contain no treated firms;
                                          // share treated differs by market arm

* Analyst's SEs treat firms as independent; compare with market-clustered SEs
reg newloan type if round==3
reg newloan type if round==3, vce(cluster survey_town)

************************************************************
* SECTION 2. CHARACTERISTICS BY TREATMENT STATUS (balance)
************************************************************

* (a) Analyst's Table 1: endline (round 3), means only
eststo clear
eststo untreated: estpost summarize part5revenue total_profit retail labor ///
    if round==3 & type==0
eststo treated:   estpost summarize part5revenue total_profit retail labor ///
    if round==3 & type==1
esttab untreated treated, cells(mean(fmt(3)) sd(par([ ]) fmt(3))) ///
    mtitle("Untreated" "Treated") label noobs

* (b) Corrected: BASELINE (round 1), regression on the four arm indicators
*     with market-clustered SEs, as in Cai-Szeidl Table 1
eststo clear
foreach x in part5revenue total_profit retail labor {
    eststo b_`x': reg `x' T50 C50 T80 C80 if round==1, vce(cluster survey_town)
}
esttab b_*, b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
    mtitle("Sales" "Profit" "Retail" "Employees") nonumbers

* Joint test: do baseline characteristics predict each arm?
local covs firmage retail labor total_profit part5revenue gender age ///
    education_college polconnection bankloan num_clients num_supplier
foreach t in T50 C50 T80 C80 {
    quietly reg `t' `covs' if round==1, vce(cluster survey_town)
    test `covs'
}

* Also: the endline gap is not a balance check (it can contain treatment effects)
reg part5revenue type if round==3, vce(cluster survey_town)
reg part5revenue type if round==1, vce(cluster survey_town)

************************************************************
* SECTION 3. DID THE TREATMENT INCREASE LOAN TAKE-UP?
************************************************************

* Peer-treatment exposure for untreated firms (as in the authors' code)
bysort survey_town round: gen townsize = _N
bysort survey_town round: egen sumtreated = total(type)
gen treatratio_peer = (sumtreated - type)/(townsize - 1)
gen inter_untpeer   = treatratio_peer*(1 - type)

* (a) Analyst: regress newloan on the market arm code (0/1/2) as if it were continuous
reg newloan survey_town_type if round==3

* (b) Corrected: firm-level treatment, market-clustered SEs
summarize newloan if round==3 & survey_town_type==0   // pure-control mean
eststo clear
eststo s3_1: reg newloan type if round==3, vce(cluster survey_town)
eststo s3_2: reg newloan type inter_untpeer if round==3, vce(cluster survey_town)
eststo s3_3: reg newloan T50 T80 C50 C80 if round==3, vce(cluster survey_town)
esttab s3_*, b(3) se(3) star(* 0.10 ** 0.05 *** 0.01)

************************************************************
* SECTION 4. TREATMENT AND NEW PRODUCT INTRODUCTION
************************************************************

gen table2_sample = round==3 & !missing(newproduct, type, retail, newloan)

* (a) Analyst's Table 2, columns 1-2 (unclustered)
reg newproduct type        if table2_sample==1
reg newproduct type retail if table2_sample==1

* Retail is uncorrelated with randomized treatment: not an omitted-variable problem
reg retail type if round==1, vce(cluster survey_town)

* (b) Corrected: clustered SEs; compare controls vs no controls
eststo clear
eststo s4_1: reg newproduct type        if table2_sample==1, vce(cluster survey_town)
eststo s4_2: reg newproduct type retail if table2_sample==1, vce(cluster survey_town)
esttab s4_1 s4_2, b(3) se(3) star(* 0.10 ** 0.05 *** 0.01)

* "Untreated" mixes pure control with untreated firms in treated markets
tab survey_town_type type if table2_sample==1, summarize(newproduct) means
eststo s4_3: reg newproduct T50 T80 C50 C80 if table2_sample==1, vce(cluster survey_town)
esttab s4_3, b(3) se(3) star(* 0.10 ** 0.05 *** 0.01)

* Authors' panel specification (direct + competitor effects, firm FE)
gen interpost  = post*type
gen inter4post = post*treatratio_comp
xtset firmid round
xtreg newproduct post interpost inter4post, fe vce(cluster survey_town)

************************************************************
* SECTION 5. DOES RECEIVING THE LOAN AFFECT PRODUCT INTRODUCTION?
************************************************************

* (a) Analyst: OLS on newloan (endogenous)
reg newproduct newloan if table2_sample==1, vce(cluster survey_town)

* (b) Corrected: randomized treatment as an instrument for newloan
*     First stage, reduced form, and 2SLS
reg newloan    type if table2_sample==1, vce(cluster survey_town)
test type                                     // first-stage F
reg newproduct type if table2_sample==1, vce(cluster survey_town)
ivregress 2sls newproduct (newloan = type) if table2_sample==1, vce(cluster survey_town)

* Robustness: add the competitor share as a second instrument
ivregress 2sls newproduct (newloan = type treatratio_comp) if table2_sample==1, ///
    vce(cluster survey_town)

log close
