
#include <stdio.h>
#include <stdlib.h>
#include "global.h"

static int processing_args=0;
static int token_ptr;
static char *line;
static int token;

char *getline();

/* This is called after an error is detected. */
scan_dump()
{
  fprintf(stderr,"scan_dump {\n");
  fprintf(stderr," processing_args = %d\n",processing_args);
  fprintf(stderr," token = %d\n",token);
  fprintf(stderr," token_ptr = %d\n",token_ptr);
  fprintf(stderr," line = \"%s\"\n",line);
  fprintf(stderr," }\n");
  }

/* This yylex switches back and forth between line by line scanning
 * and directive argument scanning.
 */
int
yylex()
{

  if (processing_args) {
#   ifdef DEBUG
    fprintf(stderr,"arg processing\n");
#   endif
    token = processargs();
    if (!token) {
      processing_args = !processing_args;
      return(T_EXPR_DONE);
      }
    }

  if (!processing_args) {
#   ifdef DEBUG
    fprintf(stderr,"line processing\n");
#   endif
    while (!(line=getline())) {
      if (!end_include()) return(0);
      }
    token = processline(line);
    switch (token) {
      case T_IF:
      case T_ELIF:
      case T_IFDEF:
      case T_IFNDEF:
        processing_args = 1;
        token_ptr = 0;
      }

  /* Output is handled here. */
    switch (token) {
      case T_FORTRAN:
        if (is_level_active()) {
          fprintf(output,"%s",line);
          }
        else {
          fprintf(output,"c_#%s",line);
          }
        break;
      default:
        fprintf(output,"c_#%s",line);
      }
    }

  /* Now that the line has been written, do include processing. */
  if (token == T_INCLUDE) {
    begin_include();
    }

# ifdef DEBUG
  print_token(token);
# endif

  return(token);
  }

int
processargs()
{
  int lenbuf=0;
  char buf[BUFSIZE];

  if ((args[token_ptr]=='\0')||(args[token_ptr]=='\n')) return(0);

  for (;;) {
    /* If I already have a bit of a token, then I may exit if I see *
     * a special character.                                         */
    if (lenbuf>0) {
      switch (args[token_ptr]) {
        case '(':
        case ')':
        case ',':
        case '=':
        case '|':
        case '&':
        case '<':
        case '>':
        case '+':
        case '-':
        case '*':
        case '/':
        case '%':
        case '!':
        case ' ':
        case '\t':
        case '\n':
        case '\0':
          buf[lenbuf] = '\0';
          if (!strcmp(buf,"defined")) return(T_DEFINED);
          if (!strcmp(buf,"streq")) return(T_STREQ);
          yylval.sval = (char *) malloc(strlen(buf)+1);
          strcpy(yylval.sval,buf);
          return(T_NAME);
          break;
        default:
          break;
        }
      }

    /* Grab a token. */
    switch(args[token_ptr]) {
      case '(':
        token_ptr++;
        return('(');
      case ')':
        token_ptr++;
        return(')');
      case ',':
        token_ptr++;
        return(',');
      case '=':
        if (args[++token_ptr]=='=') {
          token_ptr++;
          return(T_EQ);
          }
        syntax_error("Yow!  Assignment operators not allowed.  Do you mean equivalence?");
      case '|':
        if (args[++token_ptr]=='|') {
          token_ptr++;
          return(T_OR);
          }
        syntax_error("Arithmetic OR not allowed");
      case '&':
        if (args[++token_ptr]=='&') {
          token_ptr++;
          return(T_AND);
          }
        syntax_error("Arithmetic AND not allowed");
      case '<':
        if (args[++token_ptr]=='=') {
          token_ptr++;
          return(T_LESS_EQ);
          }
        return('<');
      case '>':
        if (args[++token_ptr]=='=') {
          token_ptr++;
          return(T_GREA_EQ);
          }
        return('>');
      case '+':
        token_ptr++;
        return('+');
      case '-':
        token_ptr++;
        return('-');
      case '*':
        token_ptr++;
        return('*');
      case '/':
        token_ptr++;
        return('/');
      case '%':
        token_ptr++;
        return('%');
      case '!':
        if (args[++token_ptr]=='=') {
          token_ptr++;
          return(T_NOT_EQ);
          }
        return('!');
      case ' ':
      case '\t':
        /* Skip past any initial whitespace. */
        if (lenbuf==0) token_ptr++;
        break;
      case '\n':
        if (lenbuf==0) return(0);
        break;
      default:
        buf[lenbuf++] = args[token_ptr++];
      }
    }
  }


int
processline(line)
char *line;
{
  int nwsp=0;  /* The number of whitespace chars before the directive. */
  int len;
  int token;
  char *ch;

  /* If dir and args are still around, then free them. */
  if (dir) free(dir);
  if (args) free(args);
  dir = 0;
  args = 0;

  /* None of this should be found in the input. */
  if (   !strncmp(line,"c_#",3)
      || !strncmp(line,"C_#",3)
      || !strncmp(line,"c_*BEGIN_INCLUDE",16)
      || !strncmp(line,"C_*BEGIN_INCLUDE",16)
      || !strncmp(line,"c_*BEGIN_PSI_INCLUDE",20)
      || !strncmp(line,"C_*BEGIN_PSI_INCLUDE",20)
      || !strncmp(line,"c_*END_INCLUDE",14)
      || !strncmp(line,"C_*END_INCLUDE",14)
      || !strncmp(line,"c_*BEGIN_FILE",14)
      || !strncmp(line,"C_*BEGIN_FILE",14) ) {
    syntax_error("found psipp output in an input file");
    }

  if (line[0] != '#') return(T_FORTRAN);

  /* Position ch to the first word in line following the initial '#'. */
  for (ch= &line[1]; (*ch==' ')||(*ch=='\t'); ch++) {
      nwsp++;
      }
  if (*ch=='\0') syntax_error("found a \"#\" followed by nothing\n");

  /* Compute the length of the first word. */ 
  for (len=0;
       (ch[len]!='\n')&&(ch[len]!='\0')&&(ch[len]!=' ')&&(ch[len]!='\t');
       len++);
  if (len==0) syntax_error("found a \"#\" followed by nothing\n");

  /* Copy the first word into the dir (directive) variable. */
  dir = (char *) malloc(len+1);
  strncpy(dir,ch,len);
  dir[len] = '\0';

  /* The args to the directive. */
  /* Note that all directives must be on one line, its best to not  *
   * be placing backslashes in fortran anyway.                      */
  for (ch= &line[len+1+nwsp]; (*ch==' ')||(*ch=='\t'); ch++);
  if (strlen(ch)!=0) {
    args = (char *) malloc(strlen(ch)+1);
    strcpy(args,ch);
    }
  else {
    args = NULL;
    }

  if (!strcmp(dir,"if")) {
    token = T_IF;
    }
  else if (!strcmp(dir,"elif")) {
    token = T_ELIF;
    }
  else if (!strcmp(dir,"ifdef")) {
    token = T_IFDEF;
    }
  else if (!strcmp(dir,"ifndef")) {
    token = T_IFNDEF;
    }
  else if (!strcmp(dir,"else")) {
    token = T_ELSE;
    }
  else if (!strcmp(dir,"endif")) {
    token = T_ENDIF;
    }
  else if (!strcmp(dir,"define")) {
    token = T_DEFINE;
    }
  else if (!strcmp(dir,"undef")) {
    token = T_UNDEF;
    }
  else if (!strcmp(dir,"include")) {
    token = T_INCLUDE;
    }
  else if (!strcmp(dir,"error")) {
    if (is_level_active()) {
      syntax_error("error directive encounted\n");
      }
    token = T_FORTRAN;
    }
  else syntax_error("undefined token\n");

  /* No need to process if, elif, else, and include if code is not active. */
  if ((token==T_IF)&&(!is_level_active())) token = T_INACTIVE_IF;
  if ((token==T_IFDEF)&&(!is_level_active())) token = T_INACTIVE_IF;
  if ((token==T_IFNDEF)&&(!is_level_active())) token = T_INACTIVE_IF;
  if ((token==T_ELIF)&&(!(is_previous_level_active()||is_level_active())))
    token = T_FORTRAN;
  if ((token==T_ELSE)&&(!(is_previous_level_active()||is_level_active())))
    token = T_FORTRAN;
  if ((token==T_INCLUDE)&&(!is_level_active()))
    token = T_FORTRAN;

  return token;
  }

/* Get a single line from the input, newline included.     *
 * When this routine is called, the old line is wiped out. */
char *
getline()
{
  int i, ch;
  static char buf[BUFSIZE];

  for (i=0; ; i++) {
    if (i>=BUFSIZE) {
      fprintf(stderr,"%s: input buffer size has been exceeded\n",progname);
      exit(1);
      }
    ch = getc(input);
    if (ch == EOF) {
      if (i==0) return(NULL);
      buf[i] = '\0';
      break;
      }
    else if (ch == '\n') {
      buf[i++] = ch;
      buf[i] = '\0';
      break;
      }
    buf[i] = ch;
    }
  /* Increment the global line number counter. */
  lineno++;
  return(buf);
  }
