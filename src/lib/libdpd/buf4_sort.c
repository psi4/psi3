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

  timer_on("buf4_sort");

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

      timer_on("pqsr");

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

      timer_off("pqsr");
      break;

  case prqs:

      timer_on("prqs");

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
	  

      timer_off("prqs");
      break;

  case prsq:

      timer_on("prsq");

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

      timer_off("prsq");
      break;

  case psqr:

      timer_on("psqr");

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

      timer_off("psqr");
      break;

  case psrq:

      timer_on("psrq");

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

      timer_off("psrq");
      break;

  case qprs:

      timer_on("qprs");

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

      timer_off("qprs");
      break;

  case qpsr:

      timer_on("qpsr");

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

      timer_off("qpsr");
      break;

  case qrps:
      fprintf(stderr,"\nDPD sort error: index ordering not yet coded.\n");
      dpd_error("buf_sort", stderr);
      break;

  case qrsp:
      fprintf(stderr,"\nDPD sort error: index ordering not yet coded.\n");
      dpd_error("buf_sort", stderr);
      break;

  case qspr:
      fprintf(stderr,"\nDPD sort error: index ordering not yet coded.\n");
      dpd_error("buf_sort", stderr);
      break;

  case qsrp:
      fprintf(stderr,"\nDPD sort error: index ordering not yet coded.\n");
      dpd_error("buf_sort", stderr);
      break;

  case rqps:

      timer_on("rqps");

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
      timer_off("rqps");
      break;

  case rqsp:

      timer_on("rqsp");

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

      timer_off("rqsp");
      break;

  case rpqs:

      timer_on("rpqs");

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

      timer_off("rpqs");
      break;

  case rpsq:
	  
      timer_on("rpsq");

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

      timer_off("rpsq");
      break;

  case rsqp:

      timer_on("rsqp");

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

      timer_off("rsqp");
      break;

  case rspq:

      timer_on("rspq");

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

      timer_off("rspq");
      break;

  case sqrp:

      timer_on("sqrp");

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
      timer_off("sqrp");
      break;

  case sqpr:
      fprintf(stderr,"\nDPD sort error: index ordering not yet coded.\n");
      dpd_error("buf_sort", stderr);
      break;

  case srqp:

      timer_on("srqp");

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

      timer_off("srqp");
      break;

  case srpq:
      timer_on("srpq");

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

      timer_off("srpq");
      break;

  case spqr:

      timer_on("spqr");

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
      timer_off("spqr");
      break;

  case sprq:

      timer_on("sprq");

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

      timer_off("sprq");
      break;
    }

  for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_wrt(&OutBuf, h);
      dpd_buf4_mat_irrep_close(&OutBuf, h);
      dpd_buf4_mat_irrep_close(InBuf, h);
    }

  dpd_buf4_close(&OutBuf);

  timer_off("buf4_sort");
}
