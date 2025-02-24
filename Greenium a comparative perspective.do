
encode Issuer, gen(issuer) 
drop Issuer       
encode IssuerTicker, gen(issuerticker) 
drop IssuerTicker    
destring Coupon, gen(coupon) force
drop Coupon
encode PrincipalCurrency, gen(principalcurrency) 
encode CountryofIssue, gen(countryofissue) 
encode AmountIssuedUSD, gen(amountissuedUSD)    
drop AmountIssuedUSD
encode Sector, gen(sector)   
drop Sector
encode Seniority, gen(seniority)        
drop Seniority
destring YieldtoMaturity, gen(ytm) force      
drop YieldtoMaturity
gen gbond = 0 if GreenBond == "No"
replace gbond = 1 if gbond == .
drop GreenBond
gen maturity_years = (Maturity - IssueDate)/365

drop if ytm ==.
drop if maturity_years ==.

global treatment gbond
global X coupon principalcurrency countryofissue amountissuedUSD seniority sector IssueDate maturity_years

// mean comparison with the entire dataset 
psmatch2 $treatment $X, outcome(ytm) logit
gen wmgen = _weight*ytm if _weight !=.
ttest wm, by (gbond)

//to get rid of firms that have not issued at least one green and one brown bond
bysort issuerticker: egen sumgbond = sum(gbond)
drop if sumgbond == 0
bysort issuerticker: egen countis = count(_N)
gen diff = sumgbond - countis
drop if diff == 0

//summary statistics
summarize amountissued ytm maturity coupon if gbond==1
summarize amountissued ytm maturity coupon if gbond==0
tabstat amountissued ytm maturity coupon if gbond==1, stats (median)
tabstat amountissued ytm maturity coupon if gbond==0, stats (median)
tab sector if gbond==1
tab countryofissue if gbond ==1
tab principalcurrency if gbond ==1


//mean comparison using the restricted sample
psmatch2 $treatment $X, outcome(ytm) logit
gen wm = _weight*ytm if _weight !=.

tabstat ytm, by(gbond) statistics(p5 p10 p90 p95)

//mean comparison t-test
ttest wm, by (gbond)

// Analysis for geographic area

gen america = 1 if CountryofIssue == "United States" | CountryofIssue == "Canada"
gen china = 1 if CountryofIssue == "China (Mainland)"  | CountryofIssue == "Taiwan"
gen euro = 1 if CountryofIssue == "France" | CountryofIssue == "Belgium" | CountryofIssue == "Norway" | CountryofIssue == "Sweden" | CountryofIssue == "Italy" | CountryofIssue == "Switzerland" | CountryofIssue == "Hungary" | CountryofIssue == "Iceland" | CountryofIssue == "Portugal" | CountryofIssue == "United Kingdom"
gen japan = 1 if CountryofIssue == "Japan" 
gen south_korea = 1 if CountryofIssue == "South Korea" 

replace america = 0 if america == .
replace china = 0 if china == .
replace euro = 0 if euro == .
replace japan = 0 if japan == .
replace south_korea = 0 if south_korea == .

psmatch2 $treatment $X if america == 1, outcome(ytm) logit
gen wmAM = _weight*ytm if _weight !=.

psmatch2 $treatment $X if china == 1, outcome(ytm) logit
gen wmCH = _weight*ytm if _weight !=.

psmatch2 $treatment $X if euro == 1, outcome(ytm) logit
gen wmEU = _weight*ytm if _weight !=.

psmatch2 $treatment $X if japan == 1, outcome(ytm) logit
gen wmJA = _weight*ytm if _weight !=.

psmatch2 $treatment $X if south_korea == 1, outcome(ytm) logit
gen wmSK = _weight*ytm if _weight !=.

ttest wmAM if america == 1, by (gbond) 
ttest wmCH if china == 1, by (gbond) 
ttest wmEU if euro == 1, by (gbond) 
ttest wmJA if japan == 1, by (gbond) 
ttest wmSK if south_korea == 1, by (gbond) 

// Analysis for principalcurrency

gen USD = 1 if PrincipalCurrency == "US Dollar"
gen CNH = 1 if PrincipalCurrency == "Chinese Yuan"
gen EUR = 1 if PrincipalCurrency == "Euro"
gen JPY = 1 if PrincipalCurrency == "Japanese Yen"
gen KRW = 1 if PrincipalCurrency == "South Korean Won"

replace USD = 0 if USD == .
replace CNH = 0 if CNH == .
replace EUR = 0 if EUR == .
replace JPY = 0 if JPY == .
replace KRW = 0 if KRW == .

psmatch2 $treatment $X if USD == 1, outcome(ytm) logit
gen wmUSD = _weight*ytm if _weight != .

psmatch2 $treatment $X if CNH == 1, outcome(ytm) logit
gen wmCNH = _weight*ytm if _weight !=.

psmatch2 $treatment $X if EUR == 1, outcome(ytm) logit
gen wmEUR = _weight*ytm if _weight !=.

psmatch2 $treatment $X if JPY == 1, outcome(ytm) logit
gen wmJPY = _weight*ytm if _weight !=.

psmatch2 $treatment $X if KRW == 1, outcome(ytm) logit
gen wmKRW = _weight*ytm if _weight !=.

ttest wmUSD if USD == 1, by (gbond) 
ttest wmCNH if CNH == 1, by (gbond) 
ttest wmEUR if EUR == 1, by (gbond) 
ttest wmJPY if JPY == 1, by (gbond) 
ttest wmKRW if KRW == 1, by (gbond) 

log close
