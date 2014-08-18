clear

program connect_to_mysql, plugin
plugin call connect_to_mysql, create

// Create dataset using macros passed back from plugin
local stop : word count `vars'
set obs `obs'

forvalues i = 1/`stop' {
	local var : word `i' of `vars'
	local type : word `i' of `types'
	if strpos("`type'", "str") > 0 {
		qui gen `type' `var' = ""
	}
	else {
		qui gen `type' `var' = .	
	}
}

plugin call connect_to_mysql *, load
qui compress
describe
list
set obs 4
replace id = _n in 3/4
replace var1 = 3 in 3/4
replace var2 = 4 in 3/4
replace var3 = "a" in 3/4
unab myvars: _all
foreach var of local myvars {
  local type: type `var'
  local mytypes " `mytypes' `type'"
}	
plugin call connect_to_mysql * in 3/4, write
