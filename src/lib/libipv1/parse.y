%{
#include <stdio.h>
#include <tmpl.h>
#include "ip_lib.h"
#include "ip_types.h"
#include "ip_read.gbl"
int yylex(void);

int
yywrap()
{return 1;}

int
yyerror(s)
char *s;
{ip_error(s);
 return 0;}
%}

%union {
  ip_value_t *val;
  char *str;
  }

%token <str> T_STRING
%type <val> values value array scalar
%destructor { free ($$); } T_STRING

%start input
%%

input:			group_defs
			;

group_defs:		group_defs group_def
			|
			;

group_def:		keyword ':' '(' group_defs ')'
												{ ip_pop_keyword(); }
			|	keyword ':' group_def
												{ ip_pop_keyword(); }
			|	keyword '=' value
												{ ip_assign_value($3);
												  ip_pop_keyword(); }
			;

keyword:		T_STRING
												{ ip_push_keyword($1); }
			;

value:			array
												{ $$ = $1; }
			|	scalar
												{ $$ = $1; }
			;

array:			'(' values ')' { $$ = $2; }
			;

values:			values value
												{ $$ = ip_array($1,$2); }
			|
												{ $$ = (void *)0x0; }
			;

scalar:			T_STRING
												{ $$ = ip_scalar($1); }
			;

%%

