#include <stdio.h>
#include "dpd.h"

/*
** dpd_buf_sort(): A general DPD buffer sorting function that will
** (eventually) handle all 24 possible permutations of four-index
** buffers.
**
** Arguments:
**   struct dpdbuf *InBuf: A pointer to the alread-initialized input
**     buffer.
**   int outfilenum: The PSI unit number for the target data.
**   enum indices index: The desired sorting pattern (see dpd.h).
**   int pqnum: The index combination for the bra indices for the new
**     dpd file.  N.B. this is NOT error checked for consistency
**     with the input buffer.
**   int rsnum: The index combination for the ket indices for the new
**     dpd file.  N.B. this is NOT error checked for consistency
**     with the input buffer.
**   char *label: A string labelling for this buffer.
**   int print_flag: A boolean for the print routines.
**   FILE *outfile: The formatted output file stream.
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
** -Daniel, December 1998 */

int dpd_buf_sort(struct dpdbuf *InBuf, int outfilenum, enum indices index,
		 int pqnum, int rsnum, char *label, int print_flag,
		 FILE *outfile)
{
  int h,nirreps, row, col;
  int p, q, r, s, P, Q, R, S, pq, rs, sr, pr, qs, qp, rq, qr, ps, sp, rp, sq;
  int Gp, Gq, Gr, Gs, Gpq, Grs, Gpr, Gqs, Grq, Gqr, Gps, Gsp, Grp, Gsq;
  struct dpdfile OutFile;

  nirreps = InBuf->params->nirreps;

  dpd_file_init(&OutFile, outfilenum, pqnum, rsnum, label, print_flag, outfile);

  for(h=0; h < nirreps; h++) {

      dpd_file_mat_irrep_init(&OutFile, h);

      switch(index) {
      case pqrs:
	  fprintf(outfile, "\nDPD sort error: invalid index ordering.\n");
	  dpd_error("buf_sort", outfile);
	  break;

      case pqsr:

	  /* p->p; q->q; s->r; r->s = pqsr */
      
	  dpd_buf_mat_irrep_init(InBuf, h);
	  dpd_buf_mat_irrep_rd(InBuf, h, print_flag, outfile);
      
	  for(pq=0; pq < OutFile.params->rowtot[h]; pq++) {
	      p = OutFile.params->roworb[h][pq][0];
	      q = OutFile.params->roworb[h][pq][1];

	      row = InBuf->params->rowidx[p][q];
	      
	      for(rs=0; rs < OutFile.params->coltot[h]; rs++) {
		  r = OutFile.params->colorb[h][rs][0];
		  s = OutFile.params->colorb[h][rs][1];

		  sr = InBuf->params->colidx[s][r];
	      
		  OutFile.matrix[h][pq][rs] = InBuf->matrix[h][row][sr];
		}
	    }
	  dpd_buf_mat_irrep_close(InBuf, h);
	  break;

      case prqs:

	  /* p->p; r->q; q->r; s->s = prqs */

	  for(Gp=0; Gp < nirreps; Gp++) {
	      Gq = Gp^h;
	      for(Gr=0; Gr < nirreps; Gr++) {
		  Gs = Gr^h;

		  /* Irreps on the source */
		  Gpr = Gqs = Gp^Gr;

		  dpd_buf_mat_irrep_init(InBuf, Gpr);
		  dpd_buf_mat_irrep_rd(InBuf, Gpr, print_flag, outfile);

		  for(p=0; p < OutFile.params->ppi[Gp]; p++) {
		      P = OutFile.params->poff[Gp] + p;
		      for(q=0; q < OutFile.params->qpi[Gq]; q++) {
			  Q = OutFile.params->qoff[Gq] + q;
			  pq = OutFile.params->rowidx[P][Q];

			  for(r=0; r < OutFile.params->rpi[Gr]; r++) {
			      R = OutFile.params->roff[Gr] + r;
			      pr = InBuf->params->rowidx[P][R];
			  
			      for(s=0; s < OutFile.params->spi[Gs]; s++) {
				  S = OutFile.params->soff[Gs] + s;
				  rs = OutFile.params->colidx[R][S];
				  qs = InBuf->params->colidx[Q][S];

 		  OutFile.matrix[h][pq][rs] = InBuf->matrix[Gpr][pr][qs];
			      
				}
			    }
			}
		    }

		  dpd_buf_mat_irrep_close(InBuf, Gpr);
		}
	    }
	  break;

      case prsq:

          /* p->p; r->q; s->r; q->s = psqr */

          for(Gp=0; Gp < nirreps; Gp++) {
              Gq = Gp^h;
              for(Gr=0; Gr < nirreps; Gr++) {
                  Gs = Gr^h;

                  Gps = Gqr = Gp^Gs;

                  dpd_buf_mat_irrep_init(InBuf, Gps);
                  dpd_buf_mat_irrep_rd(InBuf, Gps, print_flag, outfile);

                  for(p=0; p < OutFile.params->ppi[Gp]; p++) {
                      P = OutFile.params->poff[Gp] + p;
                      for(q=0; q < OutFile.params->qpi[Gq]; q++) {
                          Q = OutFile.params->qoff[Gq] + q;
                          pq = OutFile.params->rowidx[P][Q];

                          for(r=0; r < OutFile.params->rpi[Gr]; r++) {
                              R = OutFile.params->roff[Gr] + r;
                              qr = InBuf->params->colidx[Q][R];

                              for(s=0; s < OutFile.params->spi[Gs]; s++) {
                                  S = OutFile.params->soff[Gs] + s;
                                  rs = OutFile.params->colidx[R][S];
                                  ps = InBuf->params->rowidx[P][S];

                      OutFile.matrix[h][pq][rs] = InBuf->matrix[Gps][ps][qr];

                                }
                            }
                        }
                    }

                  dpd_buf_mat_irrep_close(InBuf, Gps);
                }
            }
	  break;

      case psqr:

	  /* p->p; s->q; q->r; r->s = prsq */
	     
	  for(Gp=0; Gp < nirreps; Gp++) {
	      Gq = Gp^h;
	      for(Gr=0; Gr < nirreps; Gr++) {
		  Gs = Gr^h;

		  Gpr = Gsq = Gp^Gr;

		  dpd_buf_mat_irrep_init(InBuf, Gpr);
		  dpd_buf_mat_irrep_rd(InBuf, Gpr, print_flag, outfile);

		  for(p=0; p < OutFile.params->ppi[Gp]; p++) {
		      P = OutFile.params->poff[Gp] + p;
		      for(q=0; q < OutFile.params->qpi[Gq]; q++) {
			  Q = OutFile.params->qoff[Gq] + q;
			  pq = OutFile.params->rowidx[P][Q];

			  for(r=0; r < OutFile.params->rpi[Gr]; r++) {
			      R = OutFile.params->roff[Gr] + r;
			      pr = InBuf->params->rowidx[P][R];

			      for(s=0; s < OutFile.params->spi[Gs]; s++) {
				  S = OutFile.params->soff[Gs] + s;
				  rs = OutFile.params->colidx[R][S];
				  sq = InBuf->params->colidx[S][Q];

		      OutFile.matrix[h][pq][rs] = InBuf->matrix[Gpr][pr][sq];
			      
				}
			    }
			}
		    }

		  dpd_buf_mat_irrep_close(InBuf, Gpr);
		}
	    }
	  break;

      case psrq:

	  /* p->p; s->q; r->r; q->s = psrq */
	     
	  for(Gp=0; Gp < nirreps; Gp++) {
	      Gq = Gp^h;
	      for(Gr=0; Gr < nirreps; Gr++) {
		  Gs = Gr^h;

		  Gps = Grq = Gp^Gs;

		  dpd_buf_mat_irrep_init(InBuf, Gps);
		  dpd_buf_mat_irrep_rd(InBuf, Gps, print_flag, outfile);

		  for(p=0; p < OutFile.params->ppi[Gp]; p++) {
		      P = OutFile.params->poff[Gp] + p;
		      for(q=0; q < OutFile.params->qpi[Gq]; q++) {
			  Q = OutFile.params->qoff[Gq] + q;
			  pq = OutFile.params->rowidx[P][Q];

			  for(r=0; r < OutFile.params->rpi[Gr]; r++) {
			      R = OutFile.params->roff[Gr] + r;
			      rq = InBuf->params->colidx[R][Q];

			      for(s=0; s < OutFile.params->spi[Gs]; s++) {
				  S = OutFile.params->soff[Gs] + s;
				  rs = OutFile.params->colidx[R][S];
				  ps = InBuf->params->rowidx[P][S];

		      OutFile.matrix[h][pq][rs] = InBuf->matrix[Gps][ps][rq];
			      
				}
			    }
			}
		    }

		  dpd_buf_mat_irrep_close(InBuf, Gps);
		}
	    }
	  break;

      case qprs:

	  /* q->p; p->q; r->r; s->s = qprs */

	  dpd_buf_mat_irrep_init(InBuf, h);
	  dpd_buf_mat_irrep_rd(InBuf, h, print_flag, outfile);

	  for(pq=0; pq < OutFile.params->rowtot[h]; pq++) {
	      p = OutFile.params->roworb[h][pq][0];
	      q = OutFile.params->roworb[h][pq][1];
	      qp = InBuf->params->rowidx[q][p];

	      for(rs=0; rs < OutFile.params->coltot[h]; rs++) {
		  r = OutFile.params->colorb[h][rs][0];
		  s = OutFile.params->colorb[h][rs][1];

		  col = InBuf->params->colidx[r][s];
		  
		  OutFile.matrix[h][pq][rs] = InBuf->matrix[h][qp][col];
		}
	    }

	  dpd_buf_mat_irrep_close(InBuf, h);
	  break;

      case qpsr:

	  /* q->p; p->q; s->r; r->s = qpsr */

	  dpd_buf_mat_irrep_init(InBuf, h);
	  dpd_buf_mat_irrep_rd(InBuf, h, print_flag, outfile);

	  for(pq=0; pq < OutFile.params->rowtot[h]; pq++) {
	      p = OutFile.params->roworb[h][pq][0];
	      q = OutFile.params->roworb[h][pq][1];
	      qp = InBuf->params->rowidx[q][p];

	      for(rs=0; rs < OutFile.params->coltot[h]; rs++) {
		  r = OutFile.params->colorb[h][rs][0];
		  s = OutFile.params->colorb[h][rs][1];
		  sr = InBuf->params->colidx[s][r];
		  
		  OutFile.matrix[h][pq][rs] = InBuf->matrix[h][qp][sr];
		}
	    }

	  dpd_buf_mat_irrep_close(InBuf, h);
	  break;

      case qrps:
	  fprintf(outfile,"\nDPD sort error: index ordering not yet coded.\n");
	  dpd_error("buf_sort", outfile);
	  break;

      case qrsp:
	  fprintf(outfile,"\nDPD sort error: index ordering not yet coded.\n");
	  dpd_error("buf_sort", outfile);
	  break;

      case qspr:
	  fprintf(outfile,"\nDPD sort error: index ordering not yet coded.\n");
	  dpd_error("buf_sort", outfile);
	  break;

      case qsrp:
	  fprintf(outfile,"\nDPD sort error: index ordering not yet coded.\n");
	  dpd_error("buf_sort", outfile);
	  break;

      case rqps:

	  /* r->p; q->q; p->r; s->s = rqps */
	  
	  for(Gp=0; Gp < nirreps; Gp++) {
	      Gq = Gp^h;
	      for(Gr=0; Gr < nirreps; Gr++) {
		  Gs = Gr^h;

		  Grq = Gps = Gp^Gs;

		  dpd_buf_mat_irrep_init(InBuf, Grq);
		  dpd_buf_mat_irrep_rd(InBuf, Grq, print_flag, outfile);

		  for(p=0; p < OutFile.params->ppi[Gp]; p++) {
		      P = OutFile.params->poff[Gp] + p;
		      for(q=0; q < OutFile.params->qpi[Gq]; q++) {
			  Q = OutFile.params->qoff[Gq] + q;
			  pq = OutFile.params->rowidx[P][Q];

			  for(r=0; r < OutFile.params->rpi[Gr]; r++) {
			      R = OutFile.params->roff[Gr] + r;
			      rq = InBuf->params->rowidx[R][Q];

			      for(s=0; s < OutFile.params->spi[Gs]; s++) {
				  S = OutFile.params->soff[Gs] + s;
				  rs = OutFile.params->colidx[R][S];
				  ps = InBuf->params->colidx[P][S];

		      OutFile.matrix[h][pq][rs] = InBuf->matrix[Grq][rq][ps];
			      
				}
			    }
			}
		    }

		  dpd_buf_mat_irrep_close(InBuf, Grq);
		}
	    }
	  break;

      case rqsp:

	  /* r->p; q->q; s->r; p->s = sqpr */
	  
	  for(Gp=0; Gp < nirreps; Gp++) {
	      Gq = Gp^h;
	      for(Gr=0; Gr < nirreps; Gr++) {
		  Gs = Gr^h;

		  Gsq = Gpr = Gp^Gr;

		  dpd_buf_mat_irrep_init(InBuf, Gsq);
		  dpd_buf_mat_irrep_rd(InBuf, Gsq, print_flag, outfile);

		  for(p=0; p < OutFile.params->ppi[Gp]; p++) {
		      P = OutFile.params->poff[Gp] + p;
		      for(q=0; q < OutFile.params->qpi[Gq]; q++) {
			  Q = OutFile.params->qoff[Gq] + q;
			  pq = OutFile.params->rowidx[P][Q];

			  for(r=0; r < OutFile.params->rpi[Gr]; r++) {
			      R = OutFile.params->roff[Gr] + r;
			      pr = InBuf->params->colidx[P][R];

			      for(s=0; s < OutFile.params->spi[Gs]; s++) {
				  S = OutFile.params->soff[Gs] + s;
				  rs = OutFile.params->colidx[R][S];
				  sq = InBuf->params->rowidx[S][Q];

		      OutFile.matrix[h][pq][rs] = InBuf->matrix[Gsq][sq][pr];
			      
				}
			    }
			}
		    }

		  dpd_buf_mat_irrep_close(InBuf, Gsq);
		}
	    }
	  break;

      case rpqs:

	  /* r->p; p->q; q->r; s->s = qrps */
	  
	  for(Gp=0; Gp < nirreps; Gp++) {
	      Gq = Gp^h;
	      for(Gr=0; Gr < nirreps; Gr++) {
		  Gs = Gr^h;

		  Gqr = Gps = Gp^Gs;

		  dpd_buf_mat_irrep_init(InBuf, Gqr);
		  dpd_buf_mat_irrep_rd(InBuf, Gqr, print_flag, outfile);

		  for(p=0; p < OutFile.params->ppi[Gp]; p++) {
		      P = OutFile.params->poff[Gp] + p;
		      for(q=0; q < OutFile.params->qpi[Gq]; q++) {
			  Q = OutFile.params->qoff[Gq] + q;
			  pq = OutFile.params->rowidx[P][Q];

			  for(r=0; r < OutFile.params->rpi[Gr]; r++) {
			      R = OutFile.params->roff[Gr] + r;
			      qr = InBuf->params->rowidx[Q][R];

			      for(s=0; s < OutFile.params->spi[Gs]; s++) {
				  S = OutFile.params->soff[Gs] + s;
				  rs = OutFile.params->colidx[R][S];
				  ps = InBuf->params->colidx[P][S];

		      OutFile.matrix[h][pq][rs] = InBuf->matrix[Gqr][qr][ps];
			      
				}
			    }
			}
		    }

		  dpd_buf_mat_irrep_close(InBuf, Gqr);
		}
	    }
	  break;

      case rpsq:
	  
	  /* r->p; p->q; s->r; q->s = qspr */
	  
	  for(Gp=0; Gp < nirreps; Gp++) {
	      Gq = Gp^h;
	      for(Gr=0; Gr < nirreps; Gr++) {
		  Gs = Gr^h;

		  Gqs = Gpr = Gp^Gr;

		  dpd_buf_mat_irrep_init(InBuf, Gqs);
		  dpd_buf_mat_irrep_rd(InBuf, Gqs, print_flag, outfile);

		  for(p=0; p < OutFile.params->ppi[Gp]; p++) {
		      P = OutFile.params->poff[Gp] + p;
		      for(q=0; q < OutFile.params->qpi[Gq]; q++) {
			  Q = OutFile.params->qoff[Gq] + q;
			  pq = OutFile.params->rowidx[P][Q];

			  for(r=0; r < OutFile.params->rpi[Gr]; r++) {
			      R = OutFile.params->roff[Gr] + r;
			      pr = InBuf->params->colidx[P][R];

			      for(s=0; s < OutFile.params->spi[Gs]; s++) {
				  S = OutFile.params->soff[Gs] + s;
				  rs = OutFile.params->colidx[R][S];
				  qs = InBuf->params->rowidx[Q][S];

		      OutFile.matrix[h][pq][rs] = InBuf->matrix[Gqs][qs][pr];
			      
				}
			    }
			}
		    }

		  dpd_buf_mat_irrep_close(InBuf, Gqs);
		}
	    }
	  break;

      case rsqp:

	  /* r->p; s->q; q->r; p->s = srpq */
	  
	  dpd_buf_mat_irrep_init(InBuf, h);
	  dpd_buf_mat_irrep_rd(InBuf, h, 0, outfile);

	  for(pq=0; pq < OutFile.params->rowtot[h]; pq++) {
	      p = OutFile.params->roworb[h][pq][0];
	      q = OutFile.params->roworb[h][pq][1];

	      col = InBuf->params->colidx[p][q];
	  
	      for(rs=0; rs < OutFile.params->coltot[h]; rs++) {
		  r = OutFile.params->colorb[h][rs][0];
		  s = OutFile.params->colorb[h][rs][1];

		  row = InBuf->params->rowidx[s][r];
		  
		  OutFile.matrix[h][pq][rs] = InBuf->matrix[h][row][col];

		}
	    }

	  dpd_buf_mat_irrep_close(InBuf, h);
	  break;

      case rspq:

	  /* r->p; s->q; p->r; q->s = rspq */
	  
	  dpd_buf_mat_irrep_init(InBuf, h);
	  dpd_buf_mat_irrep_rd(InBuf, h, 0, outfile);

	  for(pq=0; pq < OutFile.params->rowtot[h]; pq++) {
	      p = OutFile.params->roworb[h][pq][0];
	      q = OutFile.params->roworb[h][pq][1];

	      col = InBuf->params->colidx[p][q];
	  
	      for(rs=0; rs < OutFile.params->coltot[h]; rs++) {
		  r = OutFile.params->colorb[h][rs][0];
		  s = OutFile.params->colorb[h][rs][1];

		  row = InBuf->params->rowidx[r][s];
		  
		  OutFile.matrix[h][pq][rs] = InBuf->matrix[h][row][col];

		}
	    }

	  dpd_buf_mat_irrep_close(InBuf, h);
	  break;

      case sqrp:

	  /* s->p; q->q; r->r; p->s = sqrp */
	  
	  for(Gp=0; Gp < nirreps; Gp++) {
	      Gq = Gp^h;
	      for(Gr=0; Gr < nirreps; Gr++) {
		  Gs = Gr^h;
		  
		  Gsq = Grp = Gs^Gq;
		  
		  dpd_buf_mat_irrep_init(InBuf, Gsq);
		  dpd_buf_mat_irrep_rd(InBuf, Gsq, print_flag, outfile);
		  
		  for(p=0; p < OutFile.params->ppi[Gp]; p++) {
		      P = OutFile.params->poff[Gp] + p;
		      for(q=0; q < OutFile.params->qpi[Gq]; q++) {
			  Q = OutFile.params->qoff[Gq] + q;
			  pq = OutFile.params->rowidx[P][Q];
			  
			  for(r=0; r < OutFile.params->rpi[Gr]; r++) {
			      R = OutFile.params->roff[Gr] + r;
			      rp = InBuf->params->colidx[R][P];
			      
			      for(s=0; s < OutFile.params->spi[Gs]; s++) {
				  S = OutFile.params->soff[Gs] + s;
				  rs = OutFile.params->colidx[R][S];
				  sq = InBuf->params->rowidx[S][Q];
				  
		      OutFile.matrix[h][pq][rs] = InBuf->matrix[Gsq][sq][rp];
		      
				}
			    }
			}
		    }
		  dpd_buf_mat_irrep_close(InBuf, Gsq);
		}
	    }
	  break;

      case sqpr:
	  fprintf(outfile,"\nDPD sort error: index ordering not yet coded.\n");
	  dpd_error("buf_sort", outfile);
	  break;

      case srqp:

	  /* s->p; r->q; q->r; p->s = srqp */
	  
	  dpd_buf_mat_irrep_init(InBuf, h);
	  dpd_buf_mat_irrep_rd(InBuf, h, 0, outfile);

	  for(pq=0; pq < OutFile.params->rowtot[h]; pq++) {
	      p = OutFile.params->roworb[h][pq][0];
	      q = OutFile.params->roworb[h][pq][1];

	      col = InBuf->params->colidx[q][p];
	  
	      for(rs=0; rs < OutFile.params->coltot[h]; rs++) {
		  r = OutFile.params->colorb[h][rs][0];
		  s = OutFile.params->colorb[h][rs][1];

		  row = InBuf->params->rowidx[s][r];
		  
		  OutFile.matrix[h][pq][rs] = InBuf->matrix[h][row][col];

		}
	    }

	  dpd_buf_mat_irrep_close(InBuf, h);
	  break;

      case srpq:
	  fprintf(outfile,"\nDPD sort error: index ordering not yet coded.\n");
	  dpd_error("buf_sort", outfile);
	  break;

      case spqr:

	  /* s->p; p->q; q->r; r->s = qrsp */
	  
	  for(Gp=0; Gp < nirreps; Gp++) {
	      Gq = Gp^h;
	      for(Gr=0; Gr < nirreps; Gr++) {
		  Gs = Gr^h;
		  
		  Gqr = Gsp = Gq^Gr;
		  
		  dpd_buf_mat_irrep_init(InBuf, Gqr);
		  dpd_buf_mat_irrep_rd(InBuf, Gqr, print_flag, outfile);
		  
		  for(p=0; p < OutFile.params->ppi[Gp]; p++) {
		      P = OutFile.params->poff[Gp] + p;
		      for(q=0; q < OutFile.params->qpi[Gq]; q++) {
			  Q = OutFile.params->qoff[Gq] + q;
			  pq = OutFile.params->rowidx[P][Q];
			  
			  for(r=0; r < OutFile.params->rpi[Gr]; r++) {
			      R = OutFile.params->roff[Gr] + r;
			      qr = InBuf->params->rowidx[Q][R];
			      			      
			      for(s=0; s < OutFile.params->spi[Gs]; s++) {
				  S = OutFile.params->soff[Gs] + s;
				  rs = OutFile.params->colidx[R][S];
				  sp = InBuf->params->colidx[S][P];
  
		      OutFile.matrix[h][pq][rs] = InBuf->matrix[Gqr][qr][sp];
		      
				}
			    }
			}
		    }
		  dpd_buf_mat_irrep_close(InBuf, Gqr);
		}
	    }
	  break;

      case sprq:

	  /* s->p; p->q; r->r; q->s = qsrp */
	  
	  for(Gp=0; Gp < nirreps; Gp++) {
	      Gq = Gp^h;
	      for(Gr=0; Gr < nirreps; Gr++) {
		  Gs = Gr^h;
		  
		  Gqs = Grp = Gq^Gs;
		  
		  dpd_buf_mat_irrep_init(InBuf, Gqs);
		  dpd_buf_mat_irrep_rd(InBuf, Gqs, print_flag, outfile);
		  
		  for(p=0; p < OutFile.params->ppi[Gp]; p++) {
		      P = OutFile.params->poff[Gp] + p;
		      for(q=0; q < OutFile.params->qpi[Gq]; q++) {
			  Q = OutFile.params->qoff[Gq] + q;
			  pq = OutFile.params->rowidx[P][Q];
			  
			  for(r=0; r < OutFile.params->rpi[Gr]; r++) {
			      R = OutFile.params->roff[Gr] + r;
			      rp = InBuf->params->colidx[R][P];
			      
			      for(s=0; s < OutFile.params->spi[Gs]; s++) {
				  S = OutFile.params->soff[Gs] + s;
				  rs = OutFile.params->colidx[R][S];
				  qs = InBuf->params->rowidx[Q][S];
				  
		      OutFile.matrix[h][pq][rs] = InBuf->matrix[Gqs][qs][rp];
		      
				}
			    }
			}
		    }
		  dpd_buf_mat_irrep_close(InBuf, Gqs);
		}
	    }
	  break;
	}
      
      dpd_file_mat_irrep_wrt(&OutFile, h, print_flag, outfile);
      dpd_file_mat_irrep_close(&OutFile, h);
    }

  dpd_file_close(&OutFile);
}
