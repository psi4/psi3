%{
/* This has rules for parsing two different things.
 * The first is the structure of the file on a line by line basis.
 * The line by line parsing is done so that the output exactly corresponds
 * to the input, which is important when the inverse transformation to a
 * master source file is done.  The second function is to parse the
 * expressions on the if and elif directives.  It is yylex's responsibity to
 * recognize when it must feed this parser the tokens from the directive
 * statements.  This is not difficult with a hand written lexical
 * analyzer, so its really not all that complex.
 */
#include <stdio.h>
%}

%union {
  int ival;
  char *sval;
  }

%token T_FORTRAN
%token T_INCLUDE
%token T_IF
%token T_INACTIVE_IF
%token T_ELIF
%token T_IFDEF
%token T_IFNDEF
%token T_ELSE
%token T_ENDIF
%token T_DEFINE
%token T_UNDEF
%token T_EXPR_DONE
%token T_OR T_AND T_EQ T_NOT_EQ T_GREA_EQ T_LESS_EQ
%token <sval> T_NAME
%token <ival> T_DEFINED T_STREQ

%type <sval> name
%type <ival> expr defined operator

%left T_OR
%left T_AND
%left T_EQ T_NOT_EQ T_GREA_EQ T_LESS_EQ '<' '>'
%left '+' '-'
%left '*' '/' '%'
%right '!'

%start prog
%%
prog:			prog unit
			|	/* epsilon */
			;

unit:			if_unit
			|	ifdef_unit
			|	ifndef_unit
			|	inactive_if_unit
			|	define_unit
			|	undef_unit
			|	include_unit
			|	fortran
			;

include_unit:	T_INCLUDE
			;

define_unit:	T_DEFINE		{ define_name(); }
			;

undef_unit:		T_UNDEF			{ undef_name(); }
			;

ifdef_unit:		T_IFDEF			{ push_if(); }
				defined			{ if_action($3); }
				T_EXPR_DONE
				prog
				else_unit
				T_ENDIF			{ pop_if(); }
			;

ifndef_unit:	T_IFNDEF		{ push_if(); }
				defined			{ if_action(! $3); }
				T_EXPR_DONE
				prog
				else_unit
				T_ENDIF			{ pop_if(); }
			;

if_unit:		T_IF			{ push_if(); }
				expr			{ if_action($3); }
				T_EXPR_DONE
				prog
				elif_units
				else_unit
				T_ENDIF			{ pop_if(); }
			;

elif_units:		elif_units elif_unit
			|	/* epsilon */
			;

elif_unit:		T_ELIF
				expr			{ else_action($2); }
				T_EXPR_DONE
				prog
			;

else_unit:		T_ELSE			{ else_action(1); }
				prog
			|	/* epsilon */
			;

  /* Inactive if's are in code which other if's turn off. *
   * Their arguments are not looked at.                   */
inactive_if_unit:
				T_INACTIVE_IF	{ push_if(); }
				prog
				T_ENDIF			{ pop_if(); }
			;

defined:		name					{ $$ = is_name($1); }
			;

fortran:		T_FORTRAN
			;

expr:			expr '+' expr			{ $$ = $1 + $3; }
			|	expr '-' expr			{ $$ = $1 - $3; }
			|	expr '*' expr			{ $$ = $1 * $3; }
			|	expr '/' expr			{ $$ = $1 / $3; }
			|	expr '%' expr			{ $$ = $1 % $3; }
			|	expr '<' expr			{ $$ = $1 < $3; }
			|	expr '>' expr			{ $$ = $1 > $3; }
			|	expr T_LESS_EQ expr		{ $$ = $1 <= $3; }
			|	expr T_GREA_EQ expr		{ $$ = $1 >= $3; }
			|	expr T_NOT_EQ expr		{ $$ = $1 != $3; }
			|	expr T_EQ expr			{ $$ = $1 == $3; }
			|	expr T_OR expr			{ $$ = $1 || $3; }
			|	expr T_AND expr			{ $$ = $1 && $3; }
			|	'!' expr				{ $$ = !$2; }
			|	operator				{ $$ = $1; }
			|	'(' expr ')'			{ $$ = $2; }
			|	name					{ $$ = eval_name_to_int($1); }
			;

operator:		T_DEFINED '(' name ')'			{ $$ = is_name($3); }
			|	T_STREQ '(' name ',' name ')'	{ $$ = streq($3,$5); }
			;

name:			T_NAME					{ $$ = $1; }
			;

