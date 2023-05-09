/*===============================================
	Willingness-to-Pay Analysis based on Jobseeker Survey
	
===============================================*/

/*===============================================
	Setup
===============================================*/

	*=== Settings ===*
	
	clear all
	set more off 
	
	*=== Globals ===*
	
	if "`c(username)'" == "lpf0488" {
		global wif = "X:/Dropbox/Firms_Sexual_Harassment"
	}

			global analysis = "${wif}/analysis"
			global data = "${analysis}/data"
				global data_fb = "${data}/fb_survey"
				global data_isb = "${data}/student_survey"
			global output = "${analysis}/output"
			
/*===============================================
	Data Setup
===============================================*/

	*=== Import ===*

	import delimited using "${data_fb}/clean/choices_dataset.csv", clear
	
	*=== Data Setup ===*
	
	label var meta_responseid "Response Unique ID" 
		encode meta_responseid, generate(meta_responseid_n)
		xtset meta_responseid_n
	label var choiceid "Choice Number (1-10)"
	label var jobid_chosen "Letter ID of chosen job offer from that choice"
	label var jobid "Letter ID of job offer" 
	label var jobchosen "Indicator for if offer was chosen (0/1)"
		
	foreach var1 of varlist salary flexibility opportunity safety genderbalance {
		encode `var1', generate(`var1'_n)
			qui tab `var1'_n, generate(`var1'_ind)
		
		ds `var1'_ind*
			local varlist "`r(varlist)'"
			
		foreach var2 of varlist `varlist' {
			
			local substring "`var1'_n=="
			
			local label : variable label `var2'
			local substringlength : strlen local substring
			local labellength : strlen local label 
			
			local label = substr("`label'", `substringlength' + 1, `labellength')
			label var `var2' "`label'"
			
			
		}
			
	}
		// LF 5/1: Note: The ordering of these is generated automatically, so is 
		// 		not hard-coded to remain consistent. Must doublecheck with new 
		//		dataset until fixed. 
	
	*=== Macros ===* 
	
	local salary salary_ind1 salary_ind2 
	local safety safety_ind1 safety_ind3 safety_ind4 safety_ind5
	local flexibility flexibility_ind2
	local opportunity opportunity_ind2 opportunity_ind3 
	local genderbalance genderbalance_ind1 genderbalance_ind2 genderbalance_ind3 genderbalance_ind4 

	*=== Respondent Data ===*
	
	preserve 
	
		import delimited using "${data_fb}/clean/PEDL_fb_clean.csv", clear bindquote(strict) maxquotedrows(unlimited)
		keep meta_responseid a_01_age a_02_education a_03_employment b_02_gender
		
		tempfile respondents 
			save `respondents'
	
	restore 
	
	merge m:1 meta_responseid using `respondents'
		drop if _merge != 3
		
/*===============================================
	Regression
===============================================*/

// 	xtreg jobchosen salary_ind1 salary_ind2 `safety' `flexibility' `opportunity' `genderbalance', fe vce(cluster meta_responseid_n)
//	
// 	foreach var of varlist `safety' `flexibility' `opportunity' `genderbalance' {
// 		local temp = _b[`var']/(((_b[salary_ind1]/0.05) - (_b[salary_ind2]/0.05))/2)
//		
// 		display "`var' wtp: `temp'"
// 	}
//	
// 	end

	
cap program drop wtp_estimate 
program def wtp_estimate 

	syntax varlist [using], [mtitles(passthru) TITLE(passthru) replace append plain]
	
	// Full Regression
	
		qui xtreg jobchosen salary_ind1 salary_ind2 `varlist', ///
			fe vce(cluster meta_responseid_n) 
		eststo wtp
	
		foreach var of varlist `varlist' {
			
			qui estadd scalar `var'_wtp = _b[`var']/(((_b[salary_ind1]/0.05) - (_b[salary_ind2]/0.05))/2)
			
			// Test
			qui {
				nlcom _b[`var']/(((_b[salary_ind1]/0.05) - (_b[salary_ind2]/0.05))/2)
					matrix b = r(b)
					matrix V = r(V)
					local std_err = sqrt(V[1,1])
					local z = b[1,1]/`std_err'
					local pvalue = 2*normal(-abs(`z'))
			}
				
			qui estadd scalar `var'_wtp_p = `pvalue'
			
		}
	
	// By Gender
	
	foreach gen in "Male" "Female" {
		
		preserve 
		
			keep if b_02_gender == "`gen'"
			
			qui xtreg jobchosen salary_ind1 salary_ind2 `varlist', ///
				fe vce(cluster meta_responseid_n) 
			eststo wtp_`gen'
		
			foreach var of varlist `varlist' {
				
				qui estadd scalar `var'_wtp = _b[`var']/(((_b[salary_ind1]/0.05) - (_b[salary_ind2]/0.05))/2)
				
				// Test
				qui {
					nlcom _b[`var']/(((_b[salary_ind1]/0.05) - (_b[salary_ind2]/0.05))/2)
						matrix b = r(b)
						matrix V = r(V)
						local std_err = sqrt(V[1,1])
						local z = b[1,1]/`std_err'
						local pvalue = 2*normal(-abs(`z'))
				}
					
				qui estadd scalar `var'_wtp_p = `pvalue'
				
			}
		
		restore 
		
	}
	
	// Labels
	
	local stats 
	local labels
	local format
	local layout
	foreach var of varlist `varlist' {
		
		local stats `stats' `var'_wtp `var'_wtp_p
		
		local varlabel : var label `var' 		
		local labels `" `labels' "`varlabel'" "'
		local format `" `format' 3 3"'
		local layout `" `layout' "@ @" "'
		
	}
	
	local stats `stats' N N_clust
	local labels `" `labels' "# of Choices" "# of Respondents" "'
	local format `" `format' 0 0"'
	local layout `" `layout' @ "'
	
	display `" `layout' "'
	
	// Export
	
	#delimit ;
	
		esttab `using', se label `replace' `append'
			keep(salary_ind1 salary_ind2 `varlist')
			cells("b ." se(par))
			mtitles("Full Sample" "Male-Only" "Female-Only")
			stats(`stats', 
					labels(`labels') layout(`layout'))
			`title' `plain'
			note("p-values are provided in the second column" "");
			
		eststo clear;
	
	#delimit cr
	
end 

	wtp_estimate `safety' `flexibility' `opportunity' `genderbalance' using "${output}/wtp/all.csv", replace title("All") plain
	end
	
	wtp_estimate `safety' using "${output}/wtp/safety.csv", replace title("Safety Amenities")
	wtp_estimate `flexibility' using "${output}/wtp/flexibility.csv", replace title("Flexibility")
	wtp_estimate `opportunity' using "${output}/wtp/opportunity.csv", replace title("Opportunity for Promotion")
	wtp_estimate `genderbalance' using "${output}/wtp/genderbalance.csv", replace title("Gender Balance")
	wtp_estimate `safety' `flexibility' `opportunity' `genderbalance'
	
	wtp_estimate `safety' `flexibility' `opportunity' `genderbalance' using "${output}/wtp/wtp_estimates.csv", replace title("All")
	wtp_estimate `safety' using "${output}/wtp/wtp_estimates.csv", append title("Safety Amenities")
	wtp_estimate `flexibility' using "${output}/wtp/wtp_estimates.csv", append title("Flexibility")
	wtp_estimate `opportunity' using "${output}/wtp/wtp_estimates.csv", append title("Opportunity for Promotion")
	wtp_estimate `genderbalance' using "${output}/wtp/wtp_estimates.csv", append title("Gender Balance")
	
	


