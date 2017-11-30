cd ~/Desktop/Articles/JEM/src

log using logs/JEMlogs.txt, text replace name(jemlog)

import excel using ../JEMSTUDY2.xlsx, clear first case(l)

qui: replace porf = cond(porf == "P", "1", cond(porf == "F", "0", ""))
rename (inst1 studno studscor porf itemno itemjud)(inst1 stdid score pred 	 ///   
itemid itemdiff)
qui: destring pred, replace
la def passfail 0 "Fail" 1 "Pass", modify
la val pred passfail
qui: g byte tchid = cond(inst1 == "A", 1, cond(inst1 == "B", 2, 3))
drop inst1
mata:
	angoff = select(st_data(., (4..6)), rowmissing(st_data(., (4..6))) :== 0)
	cutscores = J(1, 4, .)
	for( i = 1; i <= 3; i++) {
		cutscores[1, i] = colsum(select(angoff[., 2], angoff[., 3] :== i))
	}
	cutscores[1, 4] = rowsum(cutscores[1, 1..3]) / 3
	st_matrix("cutscores", cutscores)
end

drop itemid itemdiff
forv i = 1/3 {
	g tch`i'cutscore = cutscores[1, `i']
	g tch`i'pf = score >= tch`i'cutscore
	g tch`i'dist = score - tch`i'cutscore
}
g angoff = cutscores[1, 4]
g byte angoffpf = score >= angoff
g byte angoffdist = score - angoff
foreach v of var *pf {
	la val `v' passfail
}
drop if mi(stdid)
la var stdid "Student ID #"
la var tchid "Instructor ID #"
la var score "Observed Test Score"
la var pred "Contrasting Groups Based Prediction"
la var tch1cutscore "Sum of Instructor 1 Item Difficulties"
la var tch2cutscore "Sum of Instructor 2 Item Difficulties"
la var tch3cutscore "Sum of Instructor 3 Item Difficulties"
la var angoff "Cutscore Derived through Angoff Method"

la var tch1pf "Pass/Fail Indicator Derived from Instructor 1 Total Difficulty"
la var tch2pf "Pass/Fail Indicator Derived from Instructor 1 Total Difficulty"
la var tch3pf "Pass/Fail Indicator Derived from Instructor 1 Total Difficulty"
la var angoffpf "Pass/Fail Indicator Derived using Angoff Threshold"

la var tch1dist "Difference of Observed Student Score and Instructor 1 Total Difficulty"
la var tch2dist "Difference of Observed Student Score and Instructor 1 Total Difficulty"
la var tch3dist "Difference of Observed Student Score and Instructor 1 Total Difficulty"
la var angoffdist "Difference of Observed Student Score and Angoff Threshold"

#d ;

tw 	hist score, discrete width(1) start(0) fc(yellow) lc(black) lw(vthin) ||
	kdensity score if pred == 0, lc(purple) lw(medthick) lp(solid) || 
	kdensity score if pred == 1, lc(orange) lw(medthick) lp(solid) 
	legend(region(lc(white)) nobox pos(12) rows(3) cols(1) size(small)
		label(1 "Observed Density")
		label(2 "Smoothed Predicted Failure") 
		label(3 "Smoothed Predicted Passing")) 
	yti("Density", size(small)) xti("Observed Student Score", size(small)) 
	xlab(#15, labsize(small)) xsca(range(0(1)15))
	ylab(, nogrid angle(0) labsize(small)) 
	graphr(ic(white) lc(white) fc(white)) 
	plotr(ic(white) lc(white) fc(white)) 
	ti(	"Score Distributions Based On : " 
		"Contrasting Groups Based Cutscores", c(black) size(medlarge) span)
	name(contrastingGroup, replace);
	
tw 	hist score, discrete width(1) start(0) fc(yellow) lc(black) lw(vthin) ||
 	kdensity score if angoffpf == 0, lc(purple) lw(medthick) lp(solid) || 
	kdensity score if angoffpf == 1, lc(orange) lw(medthick) lp(solid) 
	legend(region(lc(white)) nobox pos(12) rows(3) cols(1) size(small)
		label(1 "Observed Density")
		label(2 "Smoothed Predicted Failure") 
		label(3 "Smoothed Predicted Passing")) 
	yti("Density", size(small)) xti("Observed Student Score", size(small)) 
	xlab(#15, labsize(small)) xsca(range(0(1)15))
	ylab(, nogrid angle(0) labsize(small)) 
	graphr(ic(white) lc(white) fc(white)) 
	plotr(ic(white) lc(white) fc(white)) 
	ti(	"Score Distributions Based On : " 
		"Angoff Method Based Cutscores", c(black) size(medlarge) span)
	name(angoff, replace);

tw 	hist score, discrete width(1) start(0) fc(yellow) lc(black) lw(vthin) ||
 	kdensity score if tch1pf == 0, lc(purple) lw(medthick) lp(solid) || 
	kdensity score if tch1pf == 1, lc(orange) lw(medthick) lp(solid) 
	legend(region(lc(white)) nobox pos(12) rows(3) cols(1) size(small)
		label(1 "Observed Density")
		label(2 "Smoothed Predicted Failure") 
		label(3 "Smoothed Predicted Passing")) 
	yti("Density", size(small)) xti("Observed Student Score", size(small)) 
	xlab(#15, labsize(small)) xsca(range(0(1)15))
	ylab(, nogrid angle(0) labsize(small)) 
	graphr(ic(white) lc(white) fc(white)) 
	plotr(ic(white) lc(white) fc(white)) 
	ti(	"Score Distributions Based On : " 
		"Sum of Instructor 1 Estimated Item Difficulties", c(black) size(medlarge) span)
	name(tch1, replace);

tw 	hist score, discrete width(1) start(0) fc(yellow) lc(black) lw(vthin) ||
 	kdensity score if tch2pf == 0, lc(purple) lw(medthick) lp(solid) || 
	kdensity score if tch2pf == 1, lc(orange) lw(medthick) lp(solid) 
	legend(region(lc(white)) nobox pos(12) rows(3) cols(1) size(small)
		label(1 "Observed Density")
		label(2 "Smoothed Predicted Failure") 
		label(3 "Smoothed Predicted Passing")) 
	yti("Density", size(small)) xti("Observed Student Score", size(small)) 
	xlab(#15, labsize(small)) xsca(range(0(1)15))
	ylab(, nogrid angle(0) labsize(small)) 
	graphr(ic(white) lc(white) fc(white)) 
	plotr(ic(white) lc(white) fc(white)) 
	ti(	"Score Distributions Based On : " 
		"Sum of Instructor 2 Estimated Item Difficulties", c(black) size(medlarge) span)
	name(tch2, replace);

tw 	hist score, discrete width(1) start(0) fc(yellow) lc(black) lw(vthin) ||
 	kdensity score if tch3pf == 0, lc(purple) lw(medthick) lp(solid) || 
	kdensity score if tch3pf == 1, lc(orange) lw(medthick) lp(solid) 
	legend(region(lc(white)) nobox pos(12) rows(3) cols(1) size(small)
		label(1 "Observed Density")
		label(2 "Smoothed Predicted Failure") 
		label(3 "Smoothed Predicted Passing")) 
	yti("Density", size(small)) xti("Observed Student Score", size(small)) 
	xlab(#15, labsize(small)) xsca(range(0(1)15))
	ylab(, nogrid angle(0) labsize(small)) 
	graphr(ic(white) lc(white) fc(white)) 
	plotr(ic(white) lc(white) fc(white)) 
	ti(	"Score Distributions Based On : " 
		"Sum of Instructor 3 Estimated Item Difficulties", c(black) size(medlarge) span)
	name(tch3, replace);
	
gr combine contrastingGroup angoff, name(methodcombo, replace) 
	graphr(ic(white) lc(white) fc(white)) 
	plotr(ic(white) lc(white) fc(white)) ;

gr combine tch1 tch2 tch3, name(tchcombo, replace)
	graphr(ic(white) lc(white) fc(white)) 
	plotr(ic(white) lc(white) fc(white)) ;

#d cr	

foreach v in contrastingGroup angoff tch1 tch2 tch3 methodcombo tchcombo {
	qui: gr export graphs/`v'.pdf, as(pdf) name(`v') replace
}

di _n "Score Summary Statistics"
tabstat score, by(tchid) c(s) s(n mean sd min p25 p50 p75 max)

di _n "Score Distributions within Instructors by predicted outcome"

forv i = 1/3 {
	ta score pred if tchid == `i', exact
}

di _n "Score Distributions between Instructors by predicted outcome"
ta score pred, exact	


// Here is the non-parametric receiver operating curve estimation
roctab pred score, detail bamber summ

// Seems like contrasting groups yields the same result overall
// Perhaps differences between instructors

forv i = 1/3 {
	di _n(2) "Instructor ID# `i'.  Sum of Predicted Item Difficulties for this instructor = " cutscores[1, `i'] 
	roctab pred score if tchid == `i', detail bamber summ
}



log c _all
