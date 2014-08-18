// Connecting to MySQL

#include </usr/include/mysql/mysql.h>   // MySQL header file
#include "stplugin.h"			// Stata header file
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// Global MySQL results and connection	
MYSQL_RES *res ;
MYSQL mysql ;

STDLL stata_call(int argc, char *argv[])
{
	int get_results(void) ;	
	int load_data(void) ;	
	int write_data(void) ;	
	ST_retcode rc = 0;
	
	// Must specify an argument to plugin
	if (argc) {
		// Check number of arguments passed to plugin	
		if (argc > 1) {
			SF_error("too many options specified \n") ;
			return((ST_retcode) 198) ;
		}
		
		
	} else {
		SF_error("must specify the load, create, or write options \n") ;
		return((ST_retcode) 198) ;
	}

	//Option create was specified	
	if (argc == 1 && strcmp(argv[0], "create") == 0) {
		// Query database and get the data types for each variable
		rc = get_results() ;	
		if (rc) {
			mysql_free_result(res) ; 
			mysql_close(&mysql) ;
			return((ST_retcode) rc) ;
		}
	}
	//Option load was specified	
	else if (argc == 1 && strcmp(argv[0], "load") == 0) {
		// load the data from the  database
		rc = load_data() ;
		if (rc) {
			mysql_free_result(res) ; 
			mysql_close(&mysql) ;
			return((ST_retcode) rc) ;
		}
	}
	//Option write was specified	
	else if (argc == 1 && strcmp(argv[0], "write") == 0) {
		// Write data to database
		rc = write_data() ;
		if (rc) {
			return((ST_retcode) rc) ;
		}
	}
	else {
		// Invalid option 
		SF_error("option not allowed \n") ;
		return((ST_retcode) 198) ;
	}
	
	return 0;
}
	
int get_results()
{
        int num_obs, num_vars, num_bytes ;
	char query[1024] ;
	char  buff[16] ;
	char * stata_mac_vars ;
	char * stata_mac_types ;
	
	
	MYSQL_ROW row ;
	MYSQL_FIELD *field ;
		
	ST_retcode rc ;
	ST_double val = 0 ;
	// Create MySQL connection
	if(mysql_init(&mysql)==NULL)
	{
		SF_error("failed to initate MySQL connection\n") ;
		return 198 ;
	}
	// Connect to database stata_test
	if (!mysql_real_connect(&mysql,"127.0.0.1","user",
				"password","stata_test",3306,NULL,0))
	{
		SF_error("error connecting to database\n") ;
		return 198 ;
	}
	// Query database		
	strcpy(query, "SELECT * FROM stata_test") ;
	rc=mysql_real_query(&mysql,query, strlen(query)) ;
	if (rc)
	{
		SF_error("error making query\n") ;
		return 198 ;
	}
	// Get results
        res=mysql_store_result(&mysql) ;
	
	num_vars = mysql_num_fields(res) ;
	num_obs = mysql_num_rows(res) ;
	// Allocate memory for arrays
	stata_mac_vars = malloc(num_vars*33*sizeof(char)) ;
	stata_mac_types = malloc(num_vars*12) ;
	*stata_mac_vars = '\0' ;
	*stata_mac_types = '\0' ;
	// Create array of variable names and data types 	
	while((field = mysql_fetch_field(res)))
	{
		strcat(stata_mac_vars, field->name) ;	
		strcat(stata_mac_vars, " ") ;	
		if (IS_NUM(field->type)) {	// Numeric data
			strcat(stata_mac_types, "double ") ;	
		}
		else {
			strcat(stata_mac_types, "str244 ") ;	
		}
	}
	// Store variable names/types and observation number into Stata macro
	SF_macro_save("_vars", stata_mac_vars) ;
	sprintf(buff, "%i", num_obs) ;
	SF_macro_save("_obs", buff) ;
	SF_macro_save("_types", stata_mac_types) ;
	// Free memory
	free(stata_mac_vars) ;
	free(stata_mac_types) ;
	
	return 0 ;
}

int load_data()
{
	int c,r ;
	char * endp ;
	
	MYSQL_ROW row ;
	MYSQL_FIELD *field ;
	
	ST_retcode rc ;
	ST_double val = 0 ;
	// Loop over rows/columns of MySQL result set	
	for(r=1;r<=mysql_num_rows(res); r++){
		row=mysql_fetch_row(res) ;
		if(row<0) break ;
			for(c=1;c<=mysql_num_fields(res);c++){
				val = strtod(row[c-1], &endp) ;
				if (val) {	
					val = strtod(row[c-1], &endp) ;
					if(rc = SF_vstore(c, r, val)) 
					return(rc) ;
				}
				else {
					if(rc = SF_sstore(c, r, row[c-1])) 
					return(rc) ;
				}
			}
	}
	return 0 ;
}

int write_data(void)
{
	int c,r,max_len ;
	max_len = 33*SF_nvars() ;	

	char query[40+max_len] ;
	char varnames[max_len] ;
	char vartypes[max_len] ;
	char buff[16] ;
	char *var = NULL ; 	
	char *type = NULL ; 	
	char *remove_comma ;	
	
	ST_retcode rc ;
	ST_double val = 0 ;

	// Create MySQL connection
	if(mysql_init(&mysql)==NULL)
	{
		SF_error("failed to initate MySQL connection\n") ;
		return 198 ;
	}
	// Connect to database stata_test
	if (!mysql_real_connect(&mysql,"127.0.0.1","user",
				"password","stata_test",3306,NULL,0))
	{
		SF_error("error connecting to database\n") ;
		return 198 ;
	}
	// Grab variable names from Stata macro myvars
	SF_macro_use("_myvars", varnames, max_len) ;
	// Build query
	strcpy(query, "INSERT INTO stata_test (") ;
	var = strtok(varnames, " ") ;	
	// Insert varnames into query
	while (var!=NULL) {
		strcat(query, var) ;
		strcat(query, ",") ;
		var = strtok(NULL, " ") ;	
	}
	remove_comma = strrchr(query, ',');
	*remove_comma = ')';	
	strcat(query, " VALUES(") ;
	// Insert values into query
	for(r=SF_in1();r<=SF_in2(); r++){
		SF_macro_use("_mytypes", vartypes, max_len) ;
		type = strtok(vartypes, " ") ;	
		if(SF_ifobs(r)) {
			for(c=1;c<=SF_nvars();c++){
				if (strstr(type, "str")) {
					if(rc = SF_sdata(c,r,buff)) return(rc) ;
				}
				else {
					if(rc = SF_vdata(c,r,&val)) return(rc) ;
					sprintf(buff, "%g", val) ;
				}
				strcat(query, "'") ;			
				strcat(query, buff) ;			
				strcat(query, "'") ;			
				strcat(query, ",") ;
				type = strtok(NULL, " ") ;
			}
			remove_comma = strrchr(query, ',') ;
			*remove_comma = ')' ;	
			strcat(query, ",(") ;
		}
	}
	remove_comma = strrchr(query, ',');
	*remove_comma = '\0';
	// Query database
	rc=mysql_real_query(&mysql,query, strlen(query)) ;
	if (rc)
	{
		SF_error("error making query\n") ;
		return 198 ;
	}
	return 0 ;
}
