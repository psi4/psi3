/*
** TIMER.C: These functions allow one to obtain user and system
** timings for arbitrary blocks of code.  If a code block is called
** repeatedly during the course of program execution, the timer
** functions will report the block's cumulative executation time and
** the number of calls. In addition, one may time multiple code blocks
** simultaneously, and even ``overlap'' timers.
**
** To use the timer functions defined here:
**
** (1) Initialize the linked list of timers at the beginning of your
** program: timer_init();
**
** (2) Start a timer at the start of the block of code:
** timer_on("My Timer");
**
** (3) Stop the timer at the end of the block: timer_off("My Timer");
**
** (4) When all timer calls are complete, dump the linked list of
** timing data to the output file, "timer.dat": timer_done();
**
** Note that timing data is written to timer.dat only at the end of
** the timer execution, i.e., when timer_done() is called. If a code
** block is called repeatedly during the course of program execution,
** the timer functions will report the block's cumulative executation
** time and the number of calls. In addition, one may time multiple
** code blocks simultaneously, and even ``overlap'' timers (i.e., one
** timer does not need to be off when a second is started).
**
** NB this code uses system functions ctime(), time(), and times(),
** which may not quite be standard on all machines.
**
** T. Daniel Crawford, August 1999.  */

#include <stdio.h>
#include <time.h>
#include <sys/param.h>
#include <sys/times.h>

#define TIMER_KEYLEN 12
#define TIMER_OFF 0
#define TIMER_ON 1

struct timer {
    char key[TIMER_KEYLEN];
    unsigned int status;
    unsigned int calls;
    double utime;
    double stime;
    struct tms ontime;
    struct timer *next;
    struct timer *last;
};

struct timer *global_timer;
time_t timer_start, timer_end;  /* Global wall-clock on and off times */

void timer_init(void)
{
  extern struct timer *global_timer;
  extern time_t timer_start;

  timer_start = time(NULL);

  global_timer = NULL;
}

void timer_done(void)
{
  FILE *timer_out;
  extern time_t timer_start, timer_end;
  struct timer *this_timer, *next_timer;
  extern struct timer *global_timer;

  timer_end = time(NULL);

  /* Dump the timing data to timer.dat and free the timers */
  timer_out = fopen("timer.dat", "a+");
  fprintf(timer_out, "\n");
  fprintf(timer_out, "Timers On : %s", ctime(&timer_start));
  fprintf(timer_out, "Timers Off: %s", ctime(&timer_end));
  fprintf(timer_out, "\nWall Time:  %10.2f seconds\n\n",
	  (double) timer_end - timer_start);

  this_timer = global_timer;
  while(this_timer != NULL) {
      if(this_timer->calls > 1) 
	  fprintf(timer_out, "%-12s: %10.2fu %10.2fs %6d calls\n",
		  this_timer->key, this_timer->utime, this_timer->stime,
		  this_timer->calls);
      else if(this_timer->calls == 1)
	  fprintf(timer_out, "%-12s: %10.2fu %10.2fs %6d call\n",
		  this_timer->key, this_timer->utime, this_timer->stime,
		  this_timer->calls);
      next_timer = this_timer->next;
      free(this_timer);
      this_timer = next_timer;
    }

  fprintf(timer_out, 
          "\n***********************************************************\n");
  fclose(timer_out);
}

struct timer *timer_scan(char *key)
{
  extern struct timer *global_timer;
  struct timer *this_timer;

  this_timer = global_timer;
  
  while(this_timer != NULL) {
      if(!strcmp(this_timer->key,key)) return(this_timer);
      this_timer = this_timer->next;
    }

  return(this_timer);
}

struct timer *timer_last(void)
{
  extern struct timer *global_timer;
  struct timer *this_timer;

  this_timer = global_timer;

  while(this_timer != NULL) {
      if(this_timer->next == NULL) return(this_timer);
      this_timer = this_timer->next;
    }
  return(NULL);
}

void timer_on(char *key)
{
  struct timer *this_timer;

  this_timer = timer_scan(key);

  if(this_timer == NULL) { /* New timer */
      this_timer = (struct timer *) malloc(sizeof(struct timer));
      strcpy(this_timer->key,key);
      this_timer->calls = 0;
      this_timer->utime = 0;
      this_timer->stime = 0;
      this_timer->next = NULL;
      this_timer->last = timer_last();
      if(this_timer->last != NULL) this_timer->last->next = this_timer;
      else global_timer = this_timer;
    }

  if((this_timer->status == TIMER_ON) && (this_timer->calls)) {
      fprintf(stderr, "Timer %s is already on.\n", key);
      exit(1);
    }

  this_timer->status = TIMER_ON;
  this_timer->calls++;
  
  times(&(this_timer->ontime));
}

void timer_off(char *key)
{
  struct tms ontime, offtime;
  struct timer *this_timer;

  this_timer = timer_scan(key);

  if(this_timer == NULL) {
      fprintf(stderr, "Bad timer key.\n");
      exit(1);
    }

  if(this_timer->status == TIMER_OFF) {
     fprintf(stderr, "Timer %s is already off.\n", this_timer->key);
     exit(1);
    }

  ontime = this_timer->ontime;

  times(&offtime);

  this_timer->utime += ((double) (offtime.tms_utime-ontime.tms_utime))/CLK_TCK;
  this_timer->stime += ((double) (offtime.tms_stime-ontime.tms_stime))/CLK_TCK;

  this_timer->status = TIMER_OFF;
}
