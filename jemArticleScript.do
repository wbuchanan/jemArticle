// Change directory to location of source code
cd ~/Desktop/Articles/JEM/src

// Open a new log file to track any/all work done with the data 
log using logs/JEMlogs.txt, text replace name(jemlog)

// Stores division character from extended ascii range in a local macro
loc divdot `: di _char(247)'

// Store multiplication dot
loc multiply `: di _char(42)'

// Store addition operator
loc add `: di _char(43)'

// Store difference operator
loc sub `: di _char(45)'

// Store fraction/division operator
loc div `: di _char(47)'

// Loads the data set containing the predicted pass/fail indicator, observed
// scores, and item judgments
import excel using ../data/raw.xlsx, clear first case(l)

// Recodes the predicted pass/fail indicator with a numeric value
qui: replace porf = cond(porf == "P", "1", cond(porf == "F", "0", ""))

// Renames the variables from the original source
rename (inst1 studno studscor porf itemno itemjud)(inst1 stdid score pred 	 ///   
itemid itemdiff)

// Recasts the predicted pass/fail indicator as a numeric type
qui: destring pred, replace

// Defines a pass/fail value label (masking) that attributes the labels with the
// numeric values
la def passfail 0 "Fail" 1 "Pass", modify

// Applies the passfail label to the values in the pred variable/column
la val pred passfail

// Creates a numeric teacher ID (needed because Mata does not include objects 
// that store text and numeric data in a single data frame like object
qui: g byte tchid = cond(inst1 == "A", 1, cond(inst1 == "B", 2, 3))

// Starts the Mata interpreter
mata:

	// Creates the object named angoff by selecting all of the non-missing values
	// across the last three columns/variables in the data set (e.g., the item 
	// judgments from the three instructors)
	angoff = select(st_data(., (5..7)), rowmissing(st_data(., (5..7))) :== 0)
	
	// Defines a container object to store the individual instructor cutscores 
	// and the cutscore defined by averaging over these ratings
	cutscores = J(1, 4, .)
	
	// Loop over the values that identify individual instructors
	for( i = 1; i <= 3; i++) {
	
		// Replace the corresponding columns in the cutscores matrix with the 
		// sum of the item probabilities for each instructor
		cutscores[1, i] = colsum(select(angoff[., 2], angoff[., 3] :== i))
		
	} // End of Loop
	
	// Averages the sum of the item probabilities across the instructors
	cutscores[1, 4] = rowsum(cutscores[1, 1..3]) / 3
	
	// Creates a Stata matrix named cutscores that contains these values
	st_matrix("cutscores", cutscores)
	
	// Creates an object to store the ratings from rater # 1
	r1 = select(angoff[., (1..2)], angoff[., 3] :== 1)
	
	// Creates an object to store the ratings from rater # 2
	r2 = select(angoff[., (1..2)], angoff[., 3] :== 2)
	
	// Creates an object to store the ratings from rater # 3
	r3 = select(angoff[., (1..2)], angoff[., 3] :== 3)
	
	// Average item probability for instructor # 1
	mu1 = mean(r1[., 2])
	
	// Average item probability for instructor # 2
	mu2 = mean(r2[., 2])
	
	// Average item probability for instructor # 3
	mu3 = mean(r3[., 2])
	
	// Standard deviation of item probabilities for instructor # 1
	sigma1 = sqrt(variance(r1[., 2]))
	
	// Standard deviation of item probabilities for instructor # 2
	sigma2 = sqrt(variance(r2[., 2]))
	
	// Standard deviation of item probabilities for instructor # 3
	sigma3 = sqrt(variance(r3[., 2]))
	
	// Creates a container matrix to store the probabilities across columns
	items = J(3, 13, .)
	
	// Container to store the standard deviations of each item
	sigmaitems = J(1, 13, .)
	
	// Loops over the item indexes 
	for(i = 1; i <= 13; i++) {
	
		// Populates the ith column in the 1st row (Instructor 1) with the 
		// Probability of the ith item
		items[1, i] = r1[i, 2]

		// Populates the ith column in the 2nd row (Instructor 2) with the 
		// Probability of the ith item
		items[2, i] = r2[i, 2]

		// Populates the ith column in the 3rd row (Instructor 3) with the 
		// Probability of the ith item
		items[3, i] = r3[i, 2]
				
	} // End of Loop
	
	// Gets the item averages across instructors 
	muitems = mean(items)

	// Loop over items again
	for(i = 1; i <= 13; i++) {
	
		// Gets the standard deviation of item probabilities across instructors
		sigmaitems[1, i] = sqrt(variance(items[., i]))

	} // End Loop to get the standard deviations across items	
	
	// Creates the item summary matrix by appending the transpose of the average
	// and standard deviations (e.g., from row to column vectors)
	itemsummary = (muitems', sigmaitems')
	
	// Builds the rater summary matrix with average and standard deviation of 
	// each rater in the same row
	ratersummary = (mu1, sigma1 \ mu2, sigma2 \ mu3, sigma3)

	// Creates a Stata matrix named items to store the item summary matrix
	st_matrix("items", itemsummary)
	
	// Creates a Stata matrix named raters to store the rater summary matrix
	st_matrix("raters", ratersummary)
	
// End Mata interpreter	
end

// Item stem for item # 1
loc istem1 "Simplify 9`div'15"

// Item stem for item # 2
loc istem2 "Multiply and Simplify 2`div'5 `multiply' 35"

// Item stem for item # 3
loc istem3 "Multiply and Simplify 3`div'10 `multiply' 43`div'100"

// Item stem for item # 4
loc istem4 "Divide and Simplify 7`div'2 `divdot' 49`div'4"

// Item stem for item # 5
loc istem5 "Divide and Simplify 7`div'4 `divdot' 7"

// Item stem for item # 6
loc istem6 "Add and Simplify 7`div'8 `add' 7`div'8"

// Item stem for item # 7
loc istem7 "Add and Simplify 7`div'9 `add' 5`div'6"

// Item stem for item # 8
loc istem8 "Subtract and Simplify 7`div'10 `sub' 13`div'25"

// Item stem for item # 9
loc istem9 "Add (write the answer as a mixed numeral) 6 5`div'6 `add' 2 5`div'6"

// Item stem for item # 10
loc istem10 "Add 8 1`div'9 `add' 7 2`div'5"

// Item stem for item # 11
loc istem11 "Subtract 9 2`div'5 `sub' 5 1`div'3"

// Item stem for item # 12
loc istem12 "Subtract (write a mixed numeral for the answer) 27 `sub' 22 1`div'2"

// Item stem for item # 13
loc istem13 "Divide (write as mixed numeral) 12 `divdot' 1 1`div'13"

// Loops over item indices
forv i = 1/13 {

	// Contructs a macro containing the item stem and item ID 
	loc irownames `"`irownames' "Item `i' : `istem`i''""'
	//loc irownames `irownames' Item`i'

} // End of loop over item indices

// Column names for the item/rater summary matrices
loc colnms Mean "Std Dev"

// Row names for the item summary matrix
mat rownames items = `irownames'

// Row names for the rater summary matrix 
mat rownames raters = Instructor1 Instructor2 Instructor3

// Adds column names to the item summary matrix
mat colnames items = `colnms'

// Adds column names to the rater summary matrix
mat colnames raters = `colnms'

// Prints the item summary matrix to the screen
mat li items

// Prints the rater summary matrix to the screen
mat li raters

// Drops the variables containing item IDs and instructor ratings for those items
drop itemid itemdiff

// Loops over the individual teacher cutscores 
forv i = 1/3 {

	// Creates new variable containing the cut score from the ith instructor
	g tch`i'cutscore = cutscores[1, `i']
	
	// Creates a pass/fail indicator based on the cut score from the ith instructor
	g tch`i'pf = score >= tch`i'cutscore
	
	// Computes distance between observed score and cutscore from the ith instructor
	g tch`i'dist = score - tch`i'cutscore
	
} // End loop over the instructor specific cutscores

// Creates the variable angoff containing the Angoff derived cut score from all 
// instructors
g angoff = cutscores[1, 4]

// Creates a pass/fail indicator based on the Angoff derived cutscore
g byte angoffpf = score >= angoff

// Creates a variable with the distance between observed and Angoff cutscore
g byte angoffdist = score - angoff

// Drops any extra cases imported from the excel file
drop if mi(stdid)

// Gets the summary statistics of the observed score for students predicted to fail
qui: su score if pred == 0, de

// Stores the median of the scores of students predicted to fail
loc fmed `r(p50)'

// Gets the summary statistics of the observed score for students predicted to pass
qui: su score if pred == 1, de

// Stores the median of the scores of students predicted to pass
loc pmed `r(p50)'

// Computes the midpoint between the medians of the passing and failing groups
loc contgr `= ((`pmed' - `fmed') / 2) + `fmed''

// Creates a variable containing the contrasting groups cutscore
qui: g contgr = `contgr'

// Creates a pass/fail indicator based on the contrasting groups cutscore
qui: g contgrpf = score >= contgr

// Creates a variable containing the difference between the observed score and 
// the pass/fail cutscore
qui: g contgrdist = score - contgr

// Loops over all of the pass/fail indicator variables
foreach v of var *pf {

	// Applies the pass/fail value label to the variable
	la val `v' passfail
	
} // End Loop

// Adds variable labels to the variables
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
la var contgr "Contrasting Groups Cutscore"
la var contgrpf "Pass/Fail Indicator Derived from Contrasting Groups Cutscore"
la var contgrdist "Difference of Observed Student Score and Contrasting Groups Cutscore"

// Changes end of line delimiter to a semicolon
#d ;

// Creates a graph showing the distribution of the scores and smoothed densities 
// for the groups predicted to pass and fail based on instructor judgment
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
		"Instructor Predictions of Passing/Failing", c(black) size(medlarge) span)
	name(instprediction, replace);

// Creates a graph showing the distribution of the scores and smoothed densities 
// for the groups predicted to pass and fail based on Angoff pass/fail cutscore
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

// Creates a graph showing the distribution of the scores and smoothed densities 
// for the groups predicted to pass and fail based on Contrasting Groups 
// pass/fail cutscore
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
		"Contrasting Groups Based Cutscores", c(black) size(medlarge) span)
	name(contrastingGroups, replace);	

// Combines the three graphs above into a single graph	
gr combine instprediction contrastingGroups angoff, name(methodcombo, replace) 
	graphr(ic(white) lc(white) fc(white)) 
	plotr(ic(white) lc(white) fc(white)) ;

// Returns end of line delimiter to carriage return	
#d cr	

// Loops over the graphs
foreach v in instprediction contrastingGroups angoff methodcombo {

	// Exports graph to a pdf and overwrites existing version if needed
	qui: gr export graphs/`v'.pdf, as(pdf) name(`v') replace
	
} // End of Loop over graphs

// Creates a table of rater summaries
estout m(raters) using results/table1.rtf, mlabels(, none)					 ///   
ti("Table 1. Average and Standard Deviation of Raters' Ratings") replace

// Prints the item summaries to the screen
estout m(items), mlabels(, none) ti("Table 2. Average and Standard Deviation of Item Ratings") 

// Displays the summary statistics for the observed scores
estpost tabstat score, by(tchid) c(v) s(n mean sd min p25 p50 p75 max)	

// Creates a file with the results from the previous command
esttab . using results/table3.rtf, main(score) unstack wide 				 ///   
mlabels("Instructors", noti) nonum ren(count "# of Students" 				 ///   
mean "Average Score" sd "Standard Deviation" min Minimum p25 "25th %ile" 	 ///   
p50 Median p75 "75th %ile" max Maximum) noobs nonote replace 				 ///   
ti("Table 3.  Summary of Observed Scores Within and Between Instructors") 

// Prints a message to the screen
di _n "Within Instructor Comparison of Angoff and Contrasting Groups Classifications"

// Loops over the instructor IDs
forv i = 1/3 {

	// Creates two-way table based on contrasting groups and Angoff classifications
	estpost ta contgrpf angoffpf if tchid == `i', exact elabels
	
	// Exports the two-way table to a separate file
	esttab . using results/table`= 3 + `i''.rtf, cells(b) unstack 			 ///   
	varlabels(`e(labels)') eqlabels(`e(eqlabels)') collabel(, none) 		 ///   
	mlabel("Pass/Fail Derived from Angoff Cutscore") stats(p_exact, 		 ///   
	labels("Fisher Exact Test")) replace 									 ///   
	ti("Table . Contrasting Groups vs Angoff Classifications for Students in Instructor `i''s class.")
	
} // End Loop over the instructor IDs

// Prints a message to the screen
di _n "Between Instructor Comparison of Angoff and Contrasting Groups Classifications"

// Creates a two-way table showing contrasting groups vs Angoff classifications
estpost ta contgrpf angoffpf, exact elabels

// Exports the results of the two way table above 
esttab . using results/table7.rtf, cells(b) unstack varlabels(`e(labels)')   ///   
eqlabels(`e(eqlabels)') collabel(, none) replace							 ///   
mlabel("Pass/Fail Derived from Angoff Cutscore") 							 ///   
stats(p_exact, labels("Fisher Exact Test")) 								 ///   
ti("Table . Contrasting Groups vs Angoff Classifications for All Students.")

// Prints a message to the screen
di _n "Correlation of Contrasting Groups and Angoff Cut Scores and Predicted Pass/Fail"

// Estimates the tetrachoric correlation coefficients for the pass/fail indicators
// and constrains the result to be positive semi-definite
tetrachoric angoffpf contgrpf pred, posdef 

// Fits a logit model of the predicted pass/fail status on observed test scores
logit pred score

// Creates a table containing the results from a logit model fitting the predicted
// pass/fail status on the observed test score
esttab . using results/table8.rtf, label cells("b(star) se(par)") nonum		 ///   
ti("Table .  Logistic Regression of Predicted Pass/Fail on Observed Test Scores") ///   
replace

// Estimates the marginal effect of scores on pass/fail prediction probability
margins, at(score=(0(1)15)) post

// Creates a table containing the results from the estimation of the marginal 
// effects
esttab . using results/table9.rtf, cells(b(star) se(par)) rename(1._at 1 	 ///   
2._at 2 3._at 3 4._at 4 5._at 5 6._at 6 7._at 7 8._at 8 9._at 9 10._at 10 	 ///   
11._at 11 12._at 12 13._at 13 14._at 14 15._at 15 16._at 16) nonum 			 ///   
mti("Marginal Effects") collabel(, none) noobs 								 ///   
note("Standard Errors in Parentheses") replace								 ///   
ti("Table .  Predicted Probabilities by Observed Score") 

// Creates a graph of the marginal effects
marginsplot, recastci(rarea) ylab(0(0.1)1, angle(0) nogrid) 				 ///   
yti("Pr(Pass)", c(black) size(small)) graphr(ic(white) fc(white) lc(white))  ///   
plotr(ic(white) fc(white) lc(white)) ciopts(lc(black) lw(vthin) 			 ///   
fc(yellow%50)) plotop(mc(purple) mlw(vvthin) mlc(black) lc(purple) 			 ///   
lw(medium)) ti("Instructor Predicted Passing vs" "Observed Test Score", 	 ///   
c(black) size(medlarge) span) name(margplot, replace)

// Exports the graph of the marginal effects
gr export graphs/logisticMarginalEffects.pdf, as(pdf) name(margplot) replace

// Closes the log file
log c _all
