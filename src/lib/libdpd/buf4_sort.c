#include <stdio.h>
#include <stdlib.h>
#include <qt.h>
#include "dpd.h"

/*
** dpd_buf4_sort(): A general DPD buffer sorting function that will
** (eventually) handle all 24 possible permutations of four-index
** buffers.  This code assumes that all symmetry-blocks of both source and
** target buffers can be stored in core.
**
** Arguments:
**   dpdbuf4 *InBuf: A pointer to the alread-initialized input
**     buffer.
**   int outfilenum: The PSI unit number for the target data.
**   enum indices index: The desired sorting pattern (see dpd.h).
**   int pqnum: The index combination for the bra indices for the new
**     dpd file4.  N.B. this is NOT error checked for consistency
**     with the input buffer.
**   int rsnum: The index combination for the ket indices for the new
**     dpd file4.  N.B. this is NOT error checked for consistency
**     with the input buffer.
**   char *label: A string labelling for this buffer.
**
** Note that the notation for some of these multi-index sorts may be
** confusing.  When the caller requests ordering "psqr" for example,
** the desired arrangement returned is indeed Out[ps][qr] =
** In[pq][rs].  However, the index notation *within the code below*
** will appear as Out[pq][rs] = In[pr][sq]. To compute this easily,
** take the desired ordering psqr, translate its indices to p->p,
** s->q, q->r, and r->s, and finally write the source indices as
** prsq.  (This "problem" arises because I want to use the same index
** notation Out[pq][rs] for the target for every rearrangement.  For
** all pair permutations, the notation is straightforward. )
**
** Note that for pqsr, qprs, and rspq (others?), we assume that
** the InBuf and OutBuf have identical row/column orderings for all
** unswapped indices.  For example, for pqsr, we assume that the
** pq-ordering of both buffers is identical.
**
** -Daniel, December 1998
**
** Modified for new naming conventions and non-totally-symmetric data.
** TDC
** September 1999
**
** NB: Timing tests have indicated that sorts which mix bra and ket
** indices are *substantially* more expensive that those which either
** transpose bra and ket or simply mix bra or mix ket indices separately.
** The source of the problem is the multiple buf4_mat_irrep_rd() calls
** required for the b-k mixing sorts (cf. pqsr vs. prqs).  If possible, one
** should arrange contractions to use fewer b-k mixing sorts.
** TDC
** May 2000
*/

int dpd_buf4_sort(dpdbuf4 *InBuf, int outfilenum, enum indices index,
		  int pqnum, int rsnum, char *label)
{
  int h,nirreps, row, col, my_irrep, r_irrep;
  int p, q, r, s, P, Q, R, S, pq, rs, sr, pr, qs, qp, rq, qr, ps, sp, rp, sq;
  int Gp, Gq, Gr, Gs, Gpq, Grs, Gpr, Gqs, Grq, Gqr, Gps, Gsp, Grp, Gsq;
  dpdbuf4 OutBuf;

  nirreps = InBuf->params->nirreps;
  my_irrep = InBuf->file.my_irrep;

#ifdef DPD_TIMER
  timer_on("buf4_sort");
#endif

  dpd_buf4_init(&OutBuf, outfilenum, my_irrep, pqnum, rsnum,
		pqnum, rsnum, 0, label);

  /* Init input and output buffers and read in all blocks of the input */
  for(h=0; h < nirreps; h++) {

    dpd_buf4_mat_irrep_init(&OutBuf, h);

    dpd_buf4_mat_irrep_init(InBuf, h);
    dpd_buf4_mat_irrep_rd(InBuf, h);
  }


  switch(index) {
  case pqrs:
    fprintf(stderr, "\nDPD sort error: invalid index ordering.\n");
    dpd_error("buf_sort", stderr);
    break;

  case pqsr:

#ifdef DPD_TIMER
    timer_on("pqsr");
#endif

    for(h=0; h < nirreps; h++) {
      r_irrep = h^my_irrep;

      /* p->p; q->q; s->r; r->s = pqsr */
      
      for(pq=0; pq < OutBuf.params->rowtot[h]; pq++) {
	p = OutBuf.params->roworb[h][pq][0];
	q = OutBuf.params->roworb[h][pq][1];

	row = InBuf->params->rowidx[p][q];
	      
	for(rs=0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
	  r = OutBuf.params->colorb[r_irrep][rs][0];
	  s = OutBuf.params->colorb[r_irrep][rs][1];

	  sr = InBuf->params->colidx[s][r];
	      
	  OutBuf.matrix[h][pq][rs] = InBuf->matrix[h][row][sr];
	}
      }
    }

#ifdef DPD_TIMER
    timer_off("pqsr");
#endif
    break;

  case prqs:

#ifdef DPD_TIMER
    timer_on("prqs");
#endif

    /* p->p; r->q; q->r; s->s = prqs */

    for(h=0; h < nirreps; h++) {
      r_irrep = h^my_irrep;

      for(Gp=0; Gp < nirreps; Gp++) {
	Gq = Gp^h;
	for(Gr=0; Gr < nirreps; Gr++) {
	  Gs = Gr^r_irrep;

	  /* Irreps on the source */
	  Gpr = Gp^Gr;  Gqs = Gq^Gs;

	  for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
	    P = OutBuf.params->poff[Gp] + p;
	    for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
	      Q = OutBuf.params->qoff[Gq] + q;
	      pq = OutBuf.params->rowidx[P][Q];

	      for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
		R = OutBuf.params->roff[Gr] + r;
		pr = InBuf->params->rowidx[P][R];
			  
		for(s=0; s < OutBuf.params->spi[Gs]; s++) {
		  S = OutBuf.params->soff[Gs] + s;
		  rs = OutBuf.params->colidx[R][S];
		  qs = InBuf->params->colidx[Q][S];

 		  OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gpr][pr][qs];
			      
		}
	      }
	    }
	  }
	}
      }
    }
	  

#ifdef DPD_TIMER
    timer_off("prqs");
#endif
    break;

  case prsq:

#ifdef DPD_TIMER
    timer_on("prsq");
#endif

    /* p->p; r->q; s->r; q->s = psqr */

    for(h=0; h < nirreps; h++) {
      r_irrep = h^my_irrep;

      for(Gp=0; Gp < nirreps; Gp++) {
	Gq = Gp^h;
	for(Gr=0; Gr < nirreps; Gr++) {
	  Gs = Gr^r_irrep;

	  Gps = Gp^Gs;  Gqr = Gq^Gr;

	  for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
	    P = OutBuf.params->poff[Gp] + p;
	    for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
	      Q = OutBuf.params->qoff[Gq] + q;
	      pq = OutBuf.params->rowidx[P][Q];

	      for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
		R = OutBuf.params->roff[Gr] + r;
		qr = InBuf->params->colidx[Q][R];

		for(s=0; s < OutBuf.params->spi[Gs]; s++) {
		  S = OutBuf.params->soff[Gs] + s;
		  rs = OutBuf.params->colidx[R][S];
		  ps = InBuf->params->rowidx[P][S];

		  OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gps][ps][qr];

		}
	      }
	    }
	  }
	}
      }
    }

#ifdef DPD_TIMER
    timer_off("prsq");
#endif
    break;

  case psqr:

#ifdef DPD_TIMER
    timer_on("psqr");
#endif

    /* p->p; s->q; q->r; r->s = prsq */

    for(h=0; h < nirreps; h++) {
      r_irrep = h^my_irrep;
	     
      for(Gp=0; Gp < nirreps; Gp++) {
	Gq = Gp^h;
	for(Gr=0; Gr < nirreps; Gr++) {
	  Gs = Gr^r_irrep;

	  Gpr = Gp^Gr;  Gsq = Gs^Gq;

	  for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
	    P = OutBuf.params->poff[Gp] + p;
	    for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
	      Q = OutBuf.params->qoff[Gq] + q;
	      pq = OutBuf.params->rowidx[P][Q];

	      for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
		R = OutBuf.params->roff[Gr] + r;
		pr = InBuf->params->rowidx[P][R];

		for(s=0; s < OutBuf.params->spi[Gs]; s++) {
		  S = OutBuf.params->soff[Gs] + s;
		  rs = OutBuf.params->colidx[R][S];
		  sq = InBuf->params->colidx[S][Q];

		  OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gpr][pr][sq];
			      
		}
	      }
	    }
	  }
	}
      }
    }

#ifdef DPD_TIMER
    timer_off("psqr");
#endif
    break;

  case psrq:

#ifdef DPD_TIMER
    timer_on("psrq");
#endif

    /* p->p; s->q; r->r; q->s = psrq */

    for(h=0; h < nirreps; h++) {
      r_irrep = h^my_irrep;
	  
      for(Gp=0; Gp < nirreps; Gp++) {
	Gq = Gp^h;
	for(Gr=0; Gr < nirreps; Gr++) {
	  Gs = Gr^r_irrep;

	  Gps = Gp^Gs;  Grq = Gr^Gq;

	  for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
	    P = OutBuf.params->poff[Gp] + p;
	    for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
	      Q = OutBuf.params->qoff[Gq] + q;
	      pq = OutBuf.params->rowidx[P][Q];

	      for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
		R = OutBuf.params->roff[Gr] + r;
		rq = InBuf->params->colidx[R][Q];

		for(s=0; s < OutBuf.params->spi[Gs]; s++) {
		  S = OutBuf.params->soff[Gs] + s;
		  rs = OutBuf.params->colidx[R][S];
		  ps = InBuf->params->rowidx[P][S];

		  OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gps][ps][rq];
			      
		}
	      }
	    }
	  }
	}
      }
    }

#ifdef DPD_TIMER
    timer_off("psrq");
#endif
    break;

  case qprs:

#ifdef DPD_TIMER
    timer_on("qprs");
#endif

    /* q->p; p->q; r->r; s->s = qprs */

    for(h=0; h < nirreps; h++) {
      r_irrep = h^my_irrep;

      for(pq=0; pq < OutBuf.params->rowtot[h]; pq++) {
	p = OutBuf.params->roworb[h][pq][0];
	q = OutBuf.params->roworb[h][pq][1];
	qp = InBuf->params->rowidx[q][p];

	for(rs=0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
	  r = OutBuf.params->colorb[r_irrep][rs][0];
	  s = OutBuf.params->colorb[r_irrep][rs][1];

	  col = InBuf->params->colidx[r][s];
		  
	  OutBuf.matrix[h][pq][rs] = InBuf->matrix[h][qp][col];
	}
      }

    }

#ifdef DPD_TIMER
    timer_off("qprs");
#endif
    break;

  case qpsr:

#ifdef DPD_TIMER
    timer_on("qpsr");
#endif

    /* q->p; p->q; s->r; r->s = qpsr */

    for(h=0; h < nirreps; h++) {
      r_irrep = h^my_irrep;

      for(pq=0; pq < OutBuf.params->rowtot[h]; pq++) {
	p = OutBuf.params->roworb[h][pq][0];
	q = OutBuf.params->roworb[h][pq][1];
	qp = InBuf->params->rowidx[q][p];

	for(rs=0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
	  r = OutBuf.params->colorb[r_irrep][rs][0];
	  s = OutBuf.params->colorb[r_irrep][rs][1];
	  sr = InBuf->params->colidx[s][r];
		  
	  OutBuf.matrix[h][pq][rs] = InBuf->matrix[h][qp][sr];
	}
      }

    }

#ifdef DPD_TIMER
    timer_off("qpsr");
#endif
    break;

  case qrps:
    fprintf(stderr,"\nDPD sort error: qrps index ordering not yet coded.\n");
    dpd_error("buf_sort", stderr);
    break;

  case qrsp:

#ifdef DPD_TIMER
    timer_on("qrsp");
#endif

    /* q->p; r->q; s->r; p->s = spqr */

    for(h=0; h < nirreps; h++) {
      r_irrep = h^my_irrep;
	  
      for(Gp=0; Gp < nirreps; Gp++) {
	Gq = Gp^h;
	for(Gr=0; Gr < nirreps; Gr++) {
	  Gs = Gr^r_irrep;

	  Gsp = Gs^Gp; Gqr = Gq^Gr;

	  for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
	    P = OutBuf.params->poff[Gp] + p;
	    for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
	      Q = OutBuf.params->qoff[Gq] + q;
	      pq = OutBuf.params->rowidx[P][Q];

	      for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
		R = OutBuf.params->roff[Gr] + r;
		qr = InBuf->params->colidx[Q][R];

		for(s=0; s < OutBuf.params->spi[Gs]; s++) {
		  S = OutBuf.params->soff[Gs] + s;
		  rs = OutBuf.params->colidx[R][S];
		  sp = InBuf->params->rowidx[S][P];

		  OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gsp][sp][qr];
			      
		}
	      }
	    }
	  }
	}
      }
    }
#ifdef DPD_TIMER
    timer_off("qrsp");
#endif
    break;

  case qspr:
    fprintf(stderr,"\nDPD sort error: qspr index ordering not yet coded.\n");
    dpd_error("buf_sort", stderr);
    break;

  case qsrp:
    fprintf(stderr,"\nDPD sort error: qsrp index ordering not yet coded.\n");
    dpd_error("buf_sort", stderr);
    break;

  case rqps:

#ifdef DPD_TIMER
    timer_on("rqps");
#endif

    /* r->p; q->q; p->r; s->s = rqps */

    for(h=0; h < nirreps; h++) {
      r_irrep = h^my_irrep;
	  
      for(Gp=0; Gp < nirreps; Gp++) {
	Gq = Gp^h;
	for(Gr=0; Gr < nirreps; Gr++) {
	  Gs = Gr^r_irrep;

	  Grq = Gr^Gq; Gps = Gp^Gs;

	  for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
	    P = OutBuf.params->poff[Gp] + p;
	    for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
	      Q = OutBuf.params->qoff[Gq] + q;
	      pq = OutBuf.params->rowidx[P][Q];

	      for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
		R = OutBuf.params->roff[Gr] + r;
		rq = InBuf->params->rowidx[R][Q];

		for(s=0; s < OutBuf.params->spi[Gs]; s++) {
		  S = OutBuf.params->soff[Gs] + s;
		  rs = OutBuf.params->colidx[R][S];
		  ps = InBuf->params->colidx[P][S];

		  OutBuf.matrix[h][pq][rs] = InBuf->matrix[Grq][rq][ps];
			      
		}
	      }
	    }
	  }
	}
      }
    }
#ifdef DPD_TIMER
    timer_off("rqps");
#endif
    break;

  case rqsp:

#ifdef DPD_TIMER
    timer_on("rqsp");
#endif

    /* r->p; q->q; s->r; p->s = sqpr */

    for(h=0; h < nirreps; h++) {
      r_irrep = h^my_irrep;
	  
      for(Gp=0; Gp < nirreps; Gp++) {
	Gq = Gp^h;
	for(Gr=0; Gr < nirreps; Gr++) {
	  Gs = Gr^r_irrep;

	  Gsq = Gs^Gq;  Gpr = Gp^Gr;

	  for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
	    P = OutBuf.params->poff[Gp] + p;
	    for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
	      Q = OutBuf.params->qoff[Gq] + q;
	      pq = OutBuf.params->rowidx[P][Q];

	      for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
		R = OutBuf.params->roff[Gr] + r;
		pr = InBuf->params->colidx[P][R];

		for(s=0; s < OutBuf.params->spi[Gs]; s++) {
		  S = OutBuf.params->soff[Gs] + s;
		  rs = OutBuf.params->colidx[R][S];
		  sq = InBuf->params->rowidx[S][Q];

		  OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gsq][sq][pr];
			      
		}
	      }
	    }
	  }
	}
      }
    }

#ifdef DPD_TIMER
    timer_off("rqsp");
#endif
    break;

  case rpqs:

#ifdef DPD_TIMER
    timer_on("rpqs");
#endif

    /* r->p; p->q; q->r; s->s = qrps */

    for(h=0; h < nirreps; h++) {
      r_irrep = h^my_irrep;
	  
      for(Gp=0; Gp < nirreps; Gp++) {
	Gq = Gp^h;
	for(Gr=0; Gr < nirreps; Gr++) {
	  Gs = Gr^r_irrep;

	  Gqr = Gq^Gr;  Gps = Gp^Gs;

	  for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
	    P = OutBuf.params->poff[Gp] + p;
	    for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
	      Q = OutBuf.params->qoff[Gq] + q;
	      pq = OutBuf.params->rowidx[P][Q];

	      for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
		R = OutBuf.params->roff[Gr] + r;
		qr = InBuf->params->rowidx[Q][R];

		for(s=0; s < OutBuf.params->spi[Gs]; s++) {
		  S = OutBuf.params->soff[Gs] + s;
		  rs = OutBuf.params->colidx[R][S];
		  ps = InBuf->params->colidx[P][S];

		  OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gqr][qr][ps];
			      
		}
	      }
	    }
	  }
	}
      }
    }

#ifdef DPD_TIMER
    timer_off("rpqs");
#endif
    break;

  case rpsq:
	  
#ifdef DPD_TIMER
    timer_on("rpsq");
#endif

    /* r->p; p->q; s->r; q->s = qspr */

    for(h=0; h < nirreps; h++) {
      r_irrep = h^my_irrep;
	  
      for(Gp=0; Gp < nirreps; Gp++) {
	Gq = Gp^h;
	for(Gr=0; Gr < nirreps; Gr++) {
	  Gs = Gr^r_irrep;

	  Gqs = Gq^Gs;  Gpr = Gp^Gr;

	  for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
	    P = OutBuf.params->poff[Gp] + p;
	    for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
	      Q = OutBuf.params->qoff[Gq] + q;
	      pq = OutBuf.params->rowidx[P][Q];

	      for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
		R = OutBuf.params->roff[Gr] + r;
		pr = InBuf->params->colidx[P][R];

		for(s=0; s < OutBuf.params->spi[Gs]; s++) {
		  S = OutBuf.params->soff[Gs] + s;
		  rs = OutBuf.params->colidx[R][S];
		  qs = InBuf->params->rowidx[Q][S];

		  OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gqs][qs][pr];
			      
		}
	      }
	    }
	  }
	}
      }
    }

#ifdef DPD_TIMER
    timer_off("rpsq");
#endif
    break;

  case rsqp:

#ifdef DPD_TIMER
    timer_on("rsqp");
#endif

    /* r->p; s->q; q->r; p->s = srpq */

    for(h=0; h < nirreps; h++) {
      r_irrep = h^my_irrep;
	  
      for(pq=0; pq < OutBuf.params->rowtot[h]; pq++) {
	p = OutBuf.params->roworb[h][pq][0];
	q = OutBuf.params->roworb[h][pq][1];

	col = InBuf->params->colidx[p][q];
	  
	for(rs=0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
	  r = OutBuf.params->colorb[r_irrep][rs][0];
	  s = OutBuf.params->colorb[r_irrep][rs][1];

	  row = InBuf->params->rowidx[s][r];
		  
	  OutBuf.matrix[h][pq][rs] = InBuf->matrix[h][row][col];

	}
      }
    }

#ifdef DPD_TIMER
    timer_off("rsqp");
#endif
    break;

  case rspq:

#ifdef DPD_TIMER
    timer_on("rspq");
#endif

    /* r->p; s->q; p->r; q->s = rspq */

    for(h=0; h < nirreps; h++) {
      r_irrep = h^my_irrep;
	  
      for(pq=0; pq < OutBuf.params->rowtot[h]; pq++) {
	p = OutBuf.params->roworb[h][pq][0];
	q = OutBuf.params->roworb[h][pq][1];

	col = InBuf->params->colidx[p][q];
	  
	for(rs=0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
	  r = OutBuf.params->colorb[r_irrep][rs][0];
	  s = OutBuf.params->colorb[r_irrep][rs][1];

	  row = InBuf->params->rowidx[r][s];
		  
	  OutBuf.matrix[h][pq][rs] = InBuf->matrix[h][row][col];

	}
      }
    }

#ifdef DPD_TIMER
    timer_off("rspq");
#endif
    break;

  case sqrp:

#ifdef DPD_TIMER
    timer_on("sqrp");
#endif

    /* s->p; q->q; r->r; p->s = sqrp */

    for(h=0; h < nirreps; h++) {
      r_irrep = h^my_irrep;
	  
      for(Gp=0; Gp < nirreps; Gp++) {
	Gq = Gp^h;
	for(Gr=0; Gr < nirreps; Gr++) {
	  Gs = Gr^r_irrep;
		  
	  Gsq = Gs^Gq;  Grp = Gr^Gp;
		  
	  for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
	    P = OutBuf.params->poff[Gp] + p;
	    for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
	      Q = OutBuf.params->qoff[Gq] + q;
	      pq = OutBuf.params->rowidx[P][Q];
			  
	      for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
		R = OutBuf.params->roff[Gr] + r;
		rp = InBuf->params->colidx[R][P];
			      
		for(s=0; s < OutBuf.params->spi[Gs]; s++) {
		  S = OutBuf.params->soff[Gs] + s;
		  rs = OutBuf.params->colidx[R][S];
		  sq = InBuf->params->rowidx[S][Q];
				  
		  OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gsq][sq][rp];
		      
		}
	      }
	    }
	  }
	}
      }
    }
#ifdef DPD_TIMER
    timer_off("sqrp");
#endif
    break;

  case sqpr:
    fprintf(stderr,"\nDPD sort error: sqpr index ordering not yet coded.\n");
    dpd_error("buf_sort", stderr);
    break;

  case srqp:

#ifdef DPD_TIMER
    timer_on("srqp");
#endif

    /* s->p; r->q; q->r; p->s = srqp */

    for(h=0; h < nirreps; h++) {
      r_irrep = h^my_irrep;
	  
      for(pq=0; pq < OutBuf.params->rowtot[h]; pq++) {
	p = OutBuf.params->roworb[h][pq][0];
	q = OutBuf.params->roworb[h][pq][1];

	col = InBuf->params->colidx[q][p];
	  
	for(rs=0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
	  r = OutBuf.params->colorb[r_irrep][rs][0];
	  s = OutBuf.params->colorb[r_irrep][rs][1];

	  row = InBuf->params->rowidx[s][r];
		  
	  OutBuf.matrix[h][pq][rs] = InBuf->matrix[h][row][col];

	}
      }
    }

#ifdef DPD_TIMER
    timer_off("srqp");
#endif
    break;

  case srpq:
#ifdef DPD_TIMER
    timer_on("srpq");
#endif

    /* s->p; r->q; p->r; q->s = rsqp */

    for(h=0; h < nirreps; h++) {
      r_irrep = h^my_irrep;

      for(pq=0; pq < OutBuf.params->rowtot[h]; pq++) {
	p = OutBuf.params->roworb[h][pq][0];
	q = OutBuf.params->roworb[h][pq][1];

	col = InBuf->params->colidx[q][p];

	for(rs=0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
	  r = OutBuf.params->colorb[r_irrep][rs][0];
	  s = OutBuf.params->colorb[r_irrep][rs][1];

	  row = InBuf->params->rowidx[r][s];

	  OutBuf.matrix[h][pq][rs] = InBuf->matrix[h][row][col];

	}
      }
    }

#ifdef DPD_TIMER
    timer_off("srpq");
#endif
    break;

  case spqr:

#ifdef DPD_TIMER
    timer_on("spqr");
#endif

    /* s->p; p->q; q->r; r->s = qrsp */

    for(h=0; h < nirreps; h++) {
      r_irrep = h^my_irrep;
	  
      for(Gp=0; Gp < nirreps; Gp++) {
	Gq = Gp^h;
	for(Gr=0; Gr < nirreps; Gr++) {
	  Gs = Gr^r_irrep;
		  
	  Gqr = Gq^Gr;  Gsp = Gs^Gp;
		  
	  for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
	    P = OutBuf.params->poff[Gp] + p;
	    for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
	      Q = OutBuf.params->qoff[Gq] + q;
	      pq = OutBuf.params->rowidx[P][Q];
			  
	      for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
		R = OutBuf.params->roff[Gr] + r;
		qr = InBuf->params->rowidx[Q][R];
			      			      
		for(s=0; s < OutBuf.params->spi[Gs]; s++) {
		  S = OutBuf.params->soff[Gs] + s;
		  rs = OutBuf.params->colidx[R][S];
		  sp = InBuf->params->colidx[S][P];
  
		  OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gqr][qr][sp];
		      
		}
	      }
	    }
	  }
	}
      }
    }
#ifdef DPD_TIMER
    timer_off("spqr");
#endif
    break;

  case sprq:

#ifdef DPD_TIMER
    timer_on("sprq");
#endif

    /* s->p; p->q; r->r; q->s = qsrp */

    for(h=0; h < nirreps; h++) {
      r_irrep = h^my_irrep;
	  
      for(Gp=0; Gp < nirreps; Gp++) {
	Gq = Gp^h;
	for(Gr=0; Gr < nirreps; Gr++) {
	  Gs = Gr^r_irrep;
		  
	  Gqs = Gq^Gs;  Grp = Gr^Gp;
		  
	  for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
	    P = OutBuf.params->poff[Gp] + p;
	    for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
	      Q = OutBuf.params->qoff[Gq] + q;
	      pq = OutBuf.params->rowidx[P][Q];
			  
	      for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
		R = OutBuf.params->roff[Gr] + r;
		rp = InBuf->params->colidx[R][P];
			      
		for(s=0; s < OutBuf.params->spi[Gs]; s++) {
		  S = OutBuf.params->soff[Gs] + s;
		  rs = OutBuf.params->colidx[R][S];
		  qs = InBuf->params->rowidx[Q][S];
				  
		  OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gqs][qs][rp];
		      
		}
	      }
	    }
	  }
	}
      }
    }

#ifdef DPD_TIMER
    timer_off("sprq");
#endif
    break;
  }

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_wrt(&OutBuf, h);
    dpd_buf4_mat_irrep_close(&OutBuf, h);
    dpd_buf4_mat_irrep_close(InBuf, h);
  }

  dpd_buf4_close(&OutBuf);

#ifdef DPD_TIMER
  timer_off("buf4_sort");
#endif
}
