#include <stdio.h>
#include <qt.h>

enum pattern {abc, acb, cab, cba, bca, bac};

void W_sort(double ***Win, double ***Wout, int nirreps, int h, int *coltot, int **colidx, 
	    int ***colorb, int *asym, int *bsym, int *aoff, int *boff,
	    int *cpi, int *coff, int **colidx_out, enum pattern index, int sum)
{
  int Ga, Gb, Gc;
  int Gab, Gac, Gca, Gcb, Gbc, Gba;
  int A, B, C, a, b, c;
  int ab, ac, ca, cb, bc, ba;

  timer_on("W_sort");

  switch(index) {

  case abc:
    fprintf(stderr, "\nW_sort: abc pattern is invalid.\n");
    exit(2);
    break;

  case acb:
    for(Gab=0; Gab < nirreps; Gab++) {
      Gc = h ^ Gab;

      for(ab=0; ab < coltot[Gab]; ab++) {

	A = colorb[Gab][ab][0];
	B = colorb[Gab][ab][1];

	Ga = asym[A]; Gb = bsym[B];
	Gac = Ga ^ Gc;

	b = B - boff[Gb];

	for(c=0; c < cpi[Gc]; c++) {
	  C = coff[Gc] + c;

	  ac = colidx_out[A][C];

	  if(sum) Wout[Gac][ac][b] += Win[Gab][ab][c];
	  else Wout[Gac][ac][b] = Win[Gab][ab][c];
	}
      }
    }

    break;

  case cab:
    for(Gab=0; Gab < nirreps; Gab++) {
      Gc = h ^ Gab;

      for(ab=0; ab < coltot[Gab]; ab++) {

	A = colorb[Gab][ab][0];
	B = colorb[Gab][ab][1];

	Ga = asym[A]; Gb = bsym[B];
	Gca = Ga ^ Gc;

	b = B - boff[Gb];

	for(c=0; c < cpi[Gc]; c++) {
	  C = coff[Gc] + c;

	  ca = colidx_out[C][A];

	  if(sum) Wout[Gca][ca][b] += Win[Gab][ab][c];
	  else Wout[Gca][ca][b] = Win[Gab][ab][c];
	}
      }
    }

    break;

  case cba:
    for(Gab=0; Gab < nirreps; Gab++) {
      Gc = h ^ Gab;

      for(ab=0; ab < coltot[Gab]; ab++) {

	A = colorb[Gab][ab][0];
	B = colorb[Gab][ab][1];

	Ga = asym[A]; Gb = bsym[B];
	a = A - aoff[Ga];

	Gcb = Gc ^ Gb;

	for(c=0; c < cpi[Gc]; c++) {
	  C = coff[Gc] + c;

	  cb = colidx_out[C][B];

	  if(sum) Wout[Gcb][cb][a] += Win[Gab][ab][c];
	  else Wout[Gcb][cb][a] = Win[Gab][ab][c];
	}
      }
    }

    break;

  case bca:
    for(Gab=0; Gab < nirreps; Gab++) {
      Gc = h ^ Gab;

      for(ab=0; ab < coltot[Gab]; ab++) {

	A = colorb[Gab][ab][0];
	B = colorb[Gab][ab][1];

	Ga = asym[A]; Gb = bsym[B];
	a = A - aoff[Ga];

	Gbc = Gb ^ Gc;

	for(c=0; c < cpi[Gc]; c++) {
	  C = coff[Gc] + c;

	  bc = colidx_out[B][C];

	  if(sum) Wout[Gbc][bc][a] += Win[Gab][ab][c];
	  else Wout[Gbc][bc][a] = Win[Gab][ab][c];
	}
      }
    }

    break;

  case bac:
    for(Gab=0; Gab < nirreps; Gab++) {
      Gc = h ^ Gab;
      Gba = Gab;

      for(ab=0; ab < coltot[Gab]; ab++) {

	A = colorb[Gab][ab][0];
	B = colorb[Gab][ab][1];

	ba = colidx_out[B][A];

	for(c=0; c < cpi[Gc]; c++) {
	  C = coff[Gc] + c;

	  if(sum) Wout[Gba][ba][c] += Win[Gab][ab][c];
	  else Wout[Gba][ba][c] = Win[Gab][ab][c];
	}
      }
    }

    break;

  }

  timer_off("W_sort");
}
