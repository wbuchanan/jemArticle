cd ~/Desktop/Articles/JEM/src

log using logs/JEMlogs.txt, text replace name(jemlog)

import excel using ../data/raw.xlsx, clear first case(l)

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
	// Rater 1 ratings
	r1 = select(angoff[., (1..2)], angoff[., 3] :== 1)
	r2 = select(angoff[., (1..2)], angoff[., 3] :== 2)
	r3 = select(angoff[., (1..2)], angoff[., 3] :== 3)
	mu1 = mean(r1[., 2])
	mu2 = mean(r2[., 2])
	mu3 = mean(r3[., 2])
	sigma1 = sqrt(variance(r1[., 2]))
	sigma2 = sqrt(variance(r2[., 2]))
	sigma3 = sqrt(variance(r3[., 2]))
	items = J(3, 13, .)
	for(i = 1; i <= 13; i++) {
		items[1, i] = r1[i, 2]
		items[2, i] = r2[i, 2]
		items[3, i] = r3[i, 2]
	}
	muitems = mean(items)
	sigmaitems = J(1, 13, .)
	for(i = 1; i <= 13; i++) {
		sigmaitems[1, i] = sqrt(variance(items[., i]))
	}
	itemsummary = (muitems', sigmaitems')
	ratersummary = (mu1, sigma1 \ mu2, sigma2 \ mu3, sigma3)
	st_matrix("items", itemsummary)
	st_matrix("raters", ratersummary)
end

forv i = 1/13 {
	loc irownames `irownames' Item`i'
}
loc colnms Mean "Std Dev"
mat rownames items = `irownames'
mat rownames raters = Instructor1 Instructor2 Instructor3
mat colnames items = `colnms'
mat colnames raters = `colnms'
mat li items
mat li raters

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
la var pred "Contrasting Groups Prediction"
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

qui: su score if pred == 0, de
loc fmed `r(p50)'

qui: su score if pred == 1, de
loc pmed `r(p50)'

loc contgr `= ((`pmed' - `fmed') / 2) + `fmed''

qui: g contgr = `contgr'
qui: g contgrpf = score >= contgr
qui: g contgrdist = score - contgr

la val contgrpf passfail
la var contgr "Contrasting Groups Cutscore"
la var contgrpf "Pass/Fail Indicator Derived from Contrasting Groups Cutscore"
la var contgrdist "Difference of Observed Student Score and Contrasting Groups Cutscore"

#d ;

tw 	hist score, discrete width(1) start(0) fc(yellow%50) lc(black) lw(vthin) ||
	kdensity score if pred == 0, lc(purple) lw(medthick) lp(solid) || 
	kdensity score if pred == 1, lc(orange) lw(medthick) lp(solid) 
	legend(region(lc(white)) nobox pos(12) rows(3) cols(1) size(small)
		label(1 "Observed Density")
		label(2 "Smoothed Predicted Failure") 
		label(3 "Smoothed Predicted Passing")) 
	yti("Density", size(small)) xti("Observed Student Score", size(small)) 
	xlab(#15, labsize(small)) xsca(range(0(1)15))
	xline(`fmed', lc(black) lw(medthick) lp(dash))
	xline(`contgr', lc(black) lw(medthick) lp(dot))
	xline(`pmed', lc(black) lw(medthick) lp(dash))
		note(	"`fmed' = Predicted Failing Median - Dashed Black Line"
				"`contgr' = Contrasting Groups Cutscore - Dotted Black Line"
				"`pmed' = Predicted Passing Median - Dashed Black Line", c(black) size(vsmall))
	ylab(, nogrid angle(0) labsize(small)) 
	graphr(ic(white) lc(white) fc(white)) 
	plotr(ic(white) lc(white) fc(white)) 
	ti(	"Score Distributions Based On : " 
		"Contrasting Groups Based Cutscores", c(black) size(medlarge) span)
	name(contrastingGroup, replace);

tw 	hist score, discrete width(1) start(0) fc(yellow%50) lc(black) lw(vthin) ||
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

/*
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
	
gr combine tch1 tch2 tch3, name(tchcombo, replace)
	graphr(ic(white) lc(white) fc(white)) 
	plotr(ic(white) lc(white) fc(white)) ;

*/

gr combine contrastingGroup angoff, name(methodcombo, replace) 
	graphr(ic(white) lc(white) fc(white)) 
	plotr(ic(white) lc(white) fc(white)) ;
	
#d cr	

foreach v in contrastingGroup angoff methodcombo {
	qui: gr export graphs/`v'.pdf, as(pdf) name(`v') replace
}

estout m(raters) using results/table1.rtf, mlabels(, none)					 ///   
ti("Table 1. Average and Standard Deviation of Raters' Ratings") replace

estout m(items) using results/table2.rtf, mlabels(, none)					 ///   
ti("Table 2. Average and Standard Deviation of Item Ratings") replace

di _n "Score Summary Statistics"
estpost tabstat score, by(tchid) c(v) s(n mean sd min p25 p50 p75 max)

esttab . using results/table3.rtf, main(score) unstack wide mlabels("Instructors", noti) nonum 		 ///   
ren(count "# of Students" mean "Average Score" sd "Standard Deviation" 		 ///   
min Minimum p25 "25th %ile" p50 Median p75 "75th %ile" max Maximum) noobs 	 ///   
ti("Table 3.  Summary of Observed Scores Within and Between Instructors") nonote replace

di _n "Within Instructor Comparison of Angoff and Contrasting Groups Classifications"

forv i = 1/3 {
	estpost ta contgrpf angoffpf if tchid == `i', exact elabels
	
	esttab . using results/table`= 3 + `i''.rtf, cells(b) unstack varlabels(`e(labels)') eqlabels(`e(eqlabels)') collabel(, none) mlabel("Pass/Fail Derived from Angoff Cutscore") stats(p_exact, labels("Fisher Exact Test")) replace ti("Table . Contrasting Groups vs Angoff Classifications for Students in Instructor `i''s class.")
	
}

di _n "Between Instructor Comparison of Angoff and Contrasting Groups Classifications"
estpost ta contgrpf angoffpf, exact elabels

esttab . using results/table7.rtf, cells(b) unstack varlabels(`e(labels)') eqlabels(`e(eqlabels)') collabel(, none) mlabel("Pass/Fail Derived from Angoff Cutscore") stats(p_exact, labels("Fisher Exact Test")) replace ti("Table . Contrasting Groups vs Angoff Classifications for All Students.")

di _n "Correlation of Contrasting Groups and Angoff Cut Scores and Predicted Pass/Fail"
tetrachoric angoffpf contgrpf pred, posdef 

logit pred score
esttab . using results/table8.rtf, label cells("b(star) se(par)") nonum		 ///   
ti("Table .  Logistic Regression of Predicted Pass/Fail on Observed Test Scores") ///   
replace

margins, at(score=(0(1)15)) post

esttab . using results/table9.rtf, cells(b(star) se(par)) rename(1._at 1 2._at 2 3._at 3 4._at 4 5._at 5 6._at 6 7._at 7 8._at 8 9._at 9 10._at 10 11._at 11 12._at 12 13._at 13 14._at 14 15._at 15 16._at 16) nonum mti("Marginal Effects") collabel(, none) noobs note("Standard Errors in Parentheses") ti("Table .  Predicted Probabilities by Observed Score") replace

marginsplot, recastci(rarea) ylab(0(0.1)1, angle(0) nogrid) yti("Pr(Pass)", c(black) size(small)) graphr(ic(white) fc(white) lc(white)) plotr(ic(white) fc(white) lc(white)) ciopts(lc(black) lw(vthin) fc(yellow%50)) plotop(mc(purple) mlw(vvthin) mlc(black) lc(purple) lw(medium)) ti("Instructor Predicted Passing vs" "Observed Test Score", c(black) size(medlarge) span) name(margplot, replace)

gr export graphs/logisticMarginalEffects.pdf, as(pdf) name(margplot) replace

log c _all
