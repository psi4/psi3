#define EXTERN
#include "includes.h"
#include "oeprop.gbl"
#include "oeprop.h"

void compute_onecgt()
{
  int i,j,k,l,ig,jg,ibf,jbf,ib,jb,jlim,kk,ll;
  int iatom, jatom, iang, jang, i_1stbf, j_1stbf, nbfi, nbfj;
  int igmin, igmax, jgmin, jgmax;
  int ixm, iym, izm, jxm, jym, jzm, iind, jind;
  int eq_shell, eq_atoms, indmax;
  int ix, iy, iz;
  double ax, ay, az, bx, by, bz, xab, yab, zab, rab2;
  double ai, bj, gamma, over_pf, norm_pf, prod_pf, dens_pf, sdens_pf, zvec_pf;
  double px, py, pz, pax, pay, paz, pbx, pby, pbz, pcx, pcy, pcz;
  double cax, cax2, cay, cay2, caz, caz2, cbx, cbx2, cby, cby2, cbz, cbz2, rca2, rcb2;
  double cax_l, cax_l_1, cax_l_2, cay_l, cay_l_1, cay_l_2, caz_l, caz_l_1, caz_l_2;
  double cbx_l, cbx_l_1, cbx_l_2, cby_l, cby_l_1, cby_l_2, cbz_l, cbz_l_1, cbz_l_2;
  int lx1, ly1, lz1, lx2, ly2, lz2;
  double x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3;
  double tx, ty, tz, t2x, t2y, t2z;
  double x,y,z,grid_normx, grid_normy, grid_normz;
  double r,r2;
  double tmp, sum1, sum2, sum3;
  PSI_FPTR junk;
  double *Density;
  double *DX_ao, *DY_ao, *DZ_ao;
  double *DX_so, *DY_so, *DZ_so;
  double **temp_mat;
  double *XX, *YY, *ZZ;
  double **evec;
  double mxx,mxy,mxz,myy,myz,mzz;  /* Components of the second moment */
  double mxxcc,myycc,mzzcc,
         mxycc,mxzcc,myzcc;       /* Correlation corrections to
                                            quadrupole moment */
  double mxxx,myyy,mzzz,mxxy,mxxz, /* Components of the third moment */
         mxyy,myyz,mxzz,myzz,mxyz;
  double mxxxcc,myyycc,mzzzcc,mxxycc,mxxzcc,
         mxyycc,myyzcc,mxzzcc,myzzcc,mxyzcc;

  
		/* Initialize some intermediates */

  qvals = init_array(3);
  qvecs = init_matrix(3,3);
  MIX = init_box(lmax+3,lmax+3,MAXMP+1);
  MIY = init_box(lmax+3,lmax+3,MAXMP+1);
  MIZ = init_box(lmax+3,lmax+3,MAXMP+1);
  MIX[0][0][0] = MIY[0][0][0] = MIZ[0][0][0] = 1.0;
  DX_ao = init_array(natri);
  DY_ao = init_array(natri);
  DZ_ao = init_array(natri);
  if (mpmax >= 2) {
    XX = init_array(natri);
    YY = init_array(natri);
    ZZ = init_array(natri);
    MOXX = init_array(nmo);
    MOYY = init_array(nmo);
    MOZZ = init_array(nmo);
  }

  if (grid) {
    grid_pts = init_matrix(nix+1,niy+1);
    if (grid == 3) {
      grid_vecX = init_matrix(nix+1,niy+1);
      grid_vecY = init_matrix(nix+1,niy+1);
      grid_vecZ = init_matrix(nix+1,niy+1);
    }
  }
  
  /* A comment on calculation of electrostatic potentials, 
     electric field and gradient integrals - since these ints are 
     not separable into products of 3 integrals we have to pack 3 indices n,l, and 
     m (powers of x,y, and z in preexponential of the primitive) into one.
     It's done by transforming the number written as nlm in the system with base 
     lmax+1 into decimal number iind (or jind). Maximum value of these numbers is 
     lmax*(lmax+1)^2, therefore AI$ matrices are boxes 
     (lmax*(lmax+1) + 1) x (lmax*(lmax+1) + 1) x (2*lmax+1); */

  indmax = lmax*(lmax+1)*(lmax+1)+1;
  if ((grid == 1) || nuc_esp)
    AI0 = init_box(indmax,indmax,2*lmax+3);

  if (nuc_esp) {
    AIX = init_box(indmax,indmax,2*lmax+2);
    AIY = init_box(indmax,indmax,2*lmax+2);
    AIZ = init_box(indmax,indmax,2*lmax+2);
    AIXX = init_box(indmax,indmax,2*lmax+1);
    AIYY = init_box(indmax,indmax,2*lmax+1);
    AIZZ = init_box(indmax,indmax,2*lmax+1);
    AIXY = init_box(indmax,indmax,2*lmax+1);
    AIXZ = init_box(indmax,indmax,2*lmax+1);
    AIYZ = init_box(indmax,indmax,2*lmax+1);
    
    phi  = init_array(natom);
    ex   = init_array(natom);
    ey   = init_array(natom);
    ez   = init_array(natom);
    dexx = init_array(natom);
    deyy = init_array(natom);
    dezz = init_array(natom);
    dexy = init_array(natom);
    dexz = init_array(natom);
    deyz = init_array(natom);
    edens = init_array(natom);
    if (spin_prop) {
      sdens = init_array(natom);
      ahfsxx = init_array(natom);
      ahfsyy = init_array(natom);
      ahfszz = init_array(natom);
      ahfsxy = init_array(natom);
      ahfsxz = init_array(natom);
      ahfsyz = init_array(natom);
    }
  }

  dx = dy = dz = 0.0;
  dx_e = dy_e = dz_e = dx_n = dy_n = dz_n = 0.0;
  dxcc = dycc = dzcc = 0.0;
  mxx = myy = mzz = mxy = mxz = myz = 0.0;
  mxxcc = myycc = mzzcc = mxycc = mxzcc = myzcc = 0.0;
  qxx = qyy = qzz = qxy = qxz = qyz = 0.0;
  qxxcc = qyycc = qzzcc = qxycc = qxzcc = qyzcc = 0.0;
  exp_x2 = exp_y2 = exp_z2 = 0.0;
  mxxx = myyy = mzzz = mxxy = mxxz = mxyy = myyz = mxzz = myzz = mxyz = 0.0;
  mxxxcc = myyycc = mzzzcc = mxxycc = mxxzcc = 
  mxyycc = myyzcc = mxzzcc = myzzcc = mxyzcc = 0.0;
  oxxx = oyyy = ozzz = oxxy = oxxz = oxyy = oyyz = oxzz = oyzz = oxyz = 0.0;
  oxxxcc = oyyycc = ozzzcc = oxxycc = oxxzcc = 
  oxyycc = oyyzcc = oxzzcc = oyzzcc = oxyzcc = 0.0;
  darw = 0.0;
  massveloc = 0.0;

	/* I-shell loop */

  for(i=0;i<nshell;i++) {
    iatom = snuc[i] - 1;
    ax = geom[iatom][0];
    ay = geom[iatom][1];
    az = geom[iatom][2];
    iang = stype[i]-1;
    izm = 1;
    iym = iang+1;
    ixm = iym*iym;
    i_1stbf = sloc[i];
    nbfi = (iang+2)*(iang+1)/2;
    igmin = sprim[i] - 1;
    igmax = igmin + snumg[i] - 1;
    
    
    	/* J-shell loop */
    
    for(j=0;j<=i;j++) {
      jatom = snuc[j] - 1;
      bx = geom[jatom][0];
      by = geom[jatom][1];
      bz = geom[jatom][2];
      jang = stype[j]-1;
      jzm = 1;
      jym = jang+1;
      jxm = jym*jym;
      j_1stbf = sloc[j];
      nbfj = (jang+2)*(jang+1)/2;
      jgmin = sprim[j] - 1;
      jgmax = jgmin + snumg[j] - 1;
      
      eq_shell = (i == j);
      eq_atoms = (iatom == jatom);
      
      xab = ax - bx;
      yab = ay - by;
      zab = az - bz;
      rab2 = xab*xab + yab*yab + zab*zab;
      
      for(ig=igmin;ig<=igmax;ig++) {
        ai = exps[ig];
        if (eq_shell)
          jgmax = ig;
        for(jg=jgmin;jg<=jgmax;jg++) {
          bj = exps[jg];
          gamma = ai + bj;
          tmp = ai*bj*rab2/gamma;
          
          over_pf = contr[ig]*contr[jg]*exp(-tmp)*pow(M_PI/gamma,1.5);
          if (eq_shell && ig != jg)
            over_pf *= 2.0;
          
          /* Compute P and intermediates */
          
          px = (ax*ai + bx*bj)/gamma;
          py = (ay*ai + by*bj)/gamma;
          pz = (az*ai + bz*bj)/gamma;
          if (eq_atoms)
            pax = pay = paz = pbx = pby = pbz = 0.0;
          else {
            pax = px - ax;  pay = py - ay;  paz = pz - az;
            pbx = px - bx;  pby = py - by;  pbz = pz - bz;
          }

		/* Computing moment integrals */

          MI_OSrecurs(pax,pay,paz,pbx,pby,pbz,gamma,iang+2,jang+2,MAXMP);


	  /*****************************************************
	  *        Computing simple multipole moments and      *
	  *             the mass-velocity (p^4) term           *
	  *****************************************************/

	  for(ibf=0;ibf<nbfi;ibf++) {
            if (eq_shell)
              jlim = ibf + 1;
            else
              jlim = nbfj;
            for(jbf=0;jbf<jlim;jbf++) {
              ib = i_1stbf + ibf - 1; jb = j_1stbf + jbf - 1;
              dens_pf = Ptot[ioff[ib]+jb] * ((ib == jb) ? 1.0 : 2.0);
              if (fabs(dens_pf) < PFACCUTOFF) continue;

              if (corr)
                zvec_pf = (ib == jb) ? zvec[ib][ib] : 
                                       zvec[ib][jb] + zvec[jb][ib];

              lx1 = xpow_bf[iang][ibf];  lx2 = xpow_bf[jang][jbf];
	      ly1 = ypow_bf[iang][ibf];  ly2 = ypow_bf[jang][jbf];
	      lz1 = zpow_bf[iang][ibf];  lz2 = zpow_bf[jang][jbf];
	      norm_pf = norm_bf[iang][ibf]*norm_bf[jang][jbf];
	      
              x0 = MIX[lx1][lx2][0];
              y0 = MIY[ly1][ly2][0];
              z0 = MIZ[lz1][lz2][0];
              x1 = MIX[lx1][lx2][1];
              y1 = MIY[ly1][ly2][1];
              z1 = MIZ[lz1][lz2][1];

	      
	         /* One-dimensional integrals over kinetic energy operators T */
	      
	      tx = bj*(2*lx2+1)*x0 - 2*bj*bj*MIX[lx1][lx2+2][0];
	      if (lx2 >= 2)
                tx -= 0.5*lx2*(lx2-1)*MIX[lx1][lx2-2][0];
	      ty = bj*(2*ly2+1)*y0 - 2*bj*bj*MIY[ly1][ly2+2][0];
              if (ly2 >= 2)
		ty -= 0.5*ly2*(ly2-1)*MIY[ly1][ly2-2][0];
	      tz = bj*(2*lz2+1)*z0 - 2*bj*bj*MIZ[lz1][lz2+2][0];
	      if (lz2 >= 2)
		tz -= 0.5*lz2*(lz2-1)*MIZ[lz1][lz2-2][0];
	      
	         /* One-dimensional integrals over operators T^2 */

	      t2x = 4*(2*lx1+1)*(2*lx2+1)*ai*bj*x0 -
		    8*(2*lx1+1)*ai*bj*bj*MIX[lx1][lx2+2][0] -
		    8*(2*lx2+1)*ai*ai*bj*MIX[lx1+2][lx2][0] +
		    16*ai*ai*bj*bj*MIX[lx1+2][lx2+2][0];
	      if (lx1 >= 2)
		t2x += 4*lx1*(lx1-1)*bj*bj*MIX[lx1-2][lx2+2][0] -
		       2*lx1*(lx1-1)*(2*lx2+1)*bj*MIX[lx1-2][lx2][0];
	      if (lx2 >= 2)
		t2x += 4*lx2*(lx2-1)*ai*ai*MIX[lx1+2][lx2-2][0] -
		       2*(2*lx1+1)*lx2*(lx2-1)*ai*MIX[lx1][lx2-2][0];
	      if ((lx1 >= 2) && (lx2 >= 2))
		t2x += lx1*(lx1-1)*lx2*(lx2-1)*MIX[lx1-2][lx2-2][0];
	      t2y = 4*(2*ly1+1)*(2*ly2+1)*ai*bj*y0 -
		    8*(2*ly1+1)*ai*bj*bj*MIY[ly1][ly2+2][0] -
		    8*(2*ly2+1)*ai*ai*bj*MIY[ly1+2][ly2][0] +
		    16*ai*ai*bj*bj*MIY[ly1+2][ly2+2][0];
	      if (ly1 >= 2)
		t2y += 4*ly1*(ly1-1)*bj*bj*MIY[ly1-2][ly2+2][0] -
		       2*ly1*(ly1-1)*(2*ly2+1)*bj*MIY[ly1-2][ly2][0];
	      if (ly2 >= 2)
		t2y += 4*ly2*(ly2-1)*ai*ai*MIY[ly1+2][ly2-2][0] -
		       2*(2*ly1+1)*ly2*(ly2-1)*ai*MIY[ly1][ly2-2][0];
	      if ((ly1 >= 2) && (ly2 >= 2))
		t2y += ly1*(ly1-1)*ly2*(ly2-1)*MIY[ly1-2][ly2-2][0];
	      t2z = 4*(2*lz1+1)*(2*lz2+1)*ai*bj*z0 -
		    8*(2*lz1+1)*ai*bj*bj*MIZ[lz1][lz2+2][0] -
		    8*(2*lz2+1)*ai*ai*bj*MIZ[lz1+2][lz2][0] +
		    16*ai*ai*bj*bj*MIZ[lz1+2][lz2+2][0];
	      if (lz1 >= 2)
		t2z += 4*lz1*(lz1-1)*bj*bj*MIZ[lz1-2][lz2+2][0] -
		       2*lz1*(lz1-1)*(2*lz2+1)*bj*MIZ[lz1-2][lz2][0];
	      if (lz2 >= 2)
		t2z += 4*lz2*(lz2-1)*ai*ai*MIZ[lz1+2][lz2-2][0] -
		       2*(2*lz1+1)*lz2*(lz2-1)*ai*MIZ[lz1][lz2-2][0];
	      if ((lz1 >= 2) && (lz2 >= 2))
		t2z += lz1*(lz1-1)*lz2*(lz2-1)*MIZ[lz1-2][lz2-2][0];

              massveloc += norm_pf*over_pf*dens_pf*(t2x*y0*z0 + x0*t2y*z0 + x0*y0*t2z +
			   8*tx*ty*z0 + 8*tx*y0*tz + 8*x0*ty*tz);


              tmp = (x1 + px*x0) * y0 * z0 * over_pf * norm_pf;
              dx -= dens_pf * tmp;
              DX_ao[ioff[ib]+jb] -= tmp;
              if (corr)
                dxcc += zvec_pf * tmp;
              tmp = x0 * (y1 + py*y0) * z0 * over_pf * norm_pf;
              dy -= dens_pf * tmp;
              DY_ao[ioff[ib]+jb] -= tmp;
              if (corr)
                dycc += zvec_pf * tmp;
              tmp = x0 * y0 * (z1 + pz*z0) * over_pf * norm_pf;
              dz -= dens_pf * tmp;
              DZ_ao[ioff[ib]+jb] -= tmp;
              if (corr)
                dzcc += zvec_pf * tmp;
              
              if (mpmax > 1) {
                x2 = MIX[xpow_bf[iang][ibf]][xpow_bf[jang][jbf]][2];
                y2 = MIY[ypow_bf[iang][ibf]][ypow_bf[jang][jbf]][2];
                z2 = MIZ[zpow_bf[iang][ibf]][zpow_bf[jang][jbf]][2];
                
                tmp = (x2 + 2*px*x1 + px*px*x0) * y0 * z0 * over_pf * norm_pf;
                mxx -= dens_pf * tmp;
                XX[ioff[ib]+jb] += tmp;
                exp_x2 += dens_pf * tmp;
                if (corr)
                  mxxcc += zvec_pf * tmp;
                
                tmp = x0 * (y2 + 2*py*y1 + py*py*y0) * z0 * over_pf * norm_pf;
                myy -= dens_pf * tmp;
                YY[ioff[ib]+jb] += tmp;
                exp_y2 += dens_pf * tmp;
                if (corr)
                  myycc += zvec_pf * tmp;
                
                tmp = x0 * y0 * (z2 + 2*pz*z1 + pz*pz*z0) * over_pf * norm_pf;
                mzz -= dens_pf * tmp;
                ZZ[ioff[ib]+jb] += tmp;
                exp_z2 += dens_pf * tmp;
                if (corr)
                  mzzcc += zvec_pf * tmp;
                
                tmp = over_pf * z0 * (x1*y1 + x1*py + y1*px + px*py) * norm_pf;
                mxy -= dens_pf * tmp;
                if (corr)
                  mxycc += zvec_pf * tmp;
                
                tmp = over_pf * y0 * (x1*z1 + x1*pz + z1*px + px*pz) * norm_pf;
                mxz -= dens_pf * tmp;
                if (corr)
                  mxzcc += zvec_pf * tmp;
                
                tmp = over_pf * x0 * (y1*z1 + y1*pz + z1*py + py*pz) * norm_pf;
                myz -= dens_pf * tmp;
                if (corr)
                  myzcc += zvec_pf * tmp;
              }
              
              if (mpmax > 2) {
                x3 = MIX[xpow_bf[iang][ibf]][xpow_bf[jang][jbf]][3];
                y3 = MIY[ypow_bf[iang][ibf]][ypow_bf[jang][jbf]][3];
                z3 = MIZ[zpow_bf[iang][ibf]][zpow_bf[jang][jbf]][3];

                tmp = norm_pf * over_pf * y0 * z0 * 
                      (x3 + 3*px*x2 + 3*px*px*x1 + px*px*px*x0);
                mxxx -= dens_pf * tmp;
                if (corr) 
                  mxxxcc += zvec_pf * tmp;

                tmp = norm_pf * over_pf * x0 * z0 *
                      (y3 + 3*py*y2 + 3*py*py*y1 + py*py*py*y0);
                myyy -= dens_pf * tmp;
                if (corr)
                  myyycc += zvec_pf * tmp;

                tmp = norm_pf * over_pf * x0 * y0 *
                      (z3 + 3*pz*z2 + 3*pz*pz*z1 + pz*pz*pz*z0);
                mzzz -= dens_pf * tmp;
                if (corr)
                  mzzzcc += zvec_pf * tmp;
                
                tmp = norm_pf * over_pf * z0 * (y1 + py*y0) * (x2 + 2*px*x1 + px*px*x0);
                mxxy -= dens_pf * tmp;
                if (corr)
                  mxxycc += zvec_pf * tmp;
                
                tmp = norm_pf * over_pf * y0 * (z1 + pz*z0) * (x2 + 2*px*x1 + px*px*x0);
                mxxz -= dens_pf * tmp;
                if (corr)
                  mxxzcc += zvec_pf * tmp;

                tmp = norm_pf * over_pf * z0 * (x1 + px*x0) * (y2 + 2*py*y1 + py*py*y0);
                mxyy -= dens_pf * tmp;
                if (corr)
                  mxyycc += zvec_pf * tmp;
                
                tmp = norm_pf * over_pf * x0 * (z1 + pz*z0) * (y2 + 2*py*y1 + py*py*y0);
                myyz -= dens_pf * tmp;
                if (corr)
                  myyzcc += zvec_pf * tmp;

                tmp = norm_pf * over_pf * (x1 + px*x0) * y0 * (z2 + 2*pz*z1 + pz*pz*z0);
                mxzz -= dens_pf * tmp;
                if (corr)
                  mxzzcc += zvec_pf * tmp;
                
                tmp = norm_pf * over_pf * x0 * (y1 + py*y0) * (z2 + 2*pz*z1 + pz*pz*z0);
                myzz -= dens_pf * tmp;
                if (corr)
                  myzzcc += zvec_pf * tmp;
                
                tmp = norm_pf * over_pf * (x1 + px*x0) * (y1 + py*y0) * (z1 + pz*z0);
                mxyz -= dens_pf * tmp;
                if (corr)
                  mxyzcc += zvec_pf * tmp;
              }
            }
          }
	  
	  
	  
	  /********************************************************
	  *   Computing electrostatic properties (see Obara and   *
	  *   Saika paper) and electron densities at the nuclei   *
	  ********************************************************/
	  
          if (nuc_esp) {
            for(k=0;k<natom;k++) {
              pcx = px - geom[k][0];
              pcy = py - geom[k][1];
              pcz = pz - geom[k][2];
              AI_OSrecurs(pax,pay,paz,pbx,pby,pbz,pcx,pcy,pcz,gamma,iang,jang);
              for(ibf=0;ibf<nbfi;ibf++) {
                if (eq_shell)
                  jlim = ibf + 1;
                else
                  jlim = nbfj;
                for(jbf=0;jbf<jlim;jbf++) {
                  ib = i_1stbf + ibf - 1; jb = j_1stbf + jbf - 1;
		  norm_pf = norm_bf[iang][ibf]*norm_bf[jang][jbf];
		  /*---

		    ATTENTION: Normalization prefactor is included in dens_pf

		    ---*/
                  dens_pf = Ptot[ioff[ib]+jb] * ((ib == jb) ? 1.0 : 2.0) * norm_pf;
                  if (fabs(dens_pf) < PFACCUTOFF) continue;
		  if (spin_prop)
		    sdens_pf = Pspin[ioff[ib]+jb] * ((ib == jb) ? 1.0 : 2.0) * norm_pf;
                  iind = zpow_bf[iang][ibf]*izm + ypow_bf[iang][ibf]*iym +
                         xpow_bf[iang][ibf]*ixm;
                  jind = zpow_bf[jang][jbf]*jzm + ypow_bf[jang][jbf]*jym +
                         xpow_bf[jang][jbf]*jxm;

                  phi[k]  += -over_pf*AI0[iind][jind][0] * dens_pf;
                  ex[k]   += over_pf*AIX[iind][jind][0] * dens_pf;
                  ey[k]   += over_pf*AIY[iind][jind][0] * dens_pf;
                  ez[k]   += over_pf*AIZ[iind][jind][0] * dens_pf;
                  dexx[k] += over_pf*AIXX[iind][jind][0] * dens_pf;
                  deyy[k] += over_pf*AIYY[iind][jind][0] * dens_pf;
                  dezz[k] += over_pf*AIZZ[iind][jind][0] * dens_pf;
                  dexy[k] += over_pf*AIXY[iind][jind][0] * dens_pf;
                  dexz[k] += over_pf*AIXZ[iind][jind][0] * dens_pf;
                  deyz[k] += over_pf*AIYZ[iind][jind][0] * dens_pf;

		  cax = geom[k][0]-ax;  cay = geom[k][1]-ay;  caz = geom[k][2]-az;
		  rca2 = cax*cax + cay*cay + caz*caz;
		  cbx = geom[k][0]-bx;  cby = geom[k][1]-by;  cbz = geom[k][2]-bz;
		  rcb2 = cbx*cbx + cby*cby + cbz*cbz;
		  tmp = contr[ig]*contr[jg]*
			pow(cax,xpow_bf[iang][ibf])*
			pow(cay,ypow_bf[iang][ibf])*
			pow(caz,zpow_bf[iang][ibf])*
			pow(cbx,xpow_bf[jang][jbf])*
			pow(cby,ypow_bf[jang][jbf])*
			pow(cbz,zpow_bf[jang][jbf])*
			exp(-ai*rca2)*exp(-bj*rcb2);
                  if (eq_shell && ig != jg)
                    tmp *= 2.0;
		  edens[k] += tmp*dens_pf;
		  if (spin_prop) {
		    sdens[k] += tmp*sdens_pf;
		    ahfsxx[k] += over_pf*AIXX[iind][jind][0] * sdens_pf;
                    ahfsyy[k] += over_pf*AIYY[iind][jind][0] * sdens_pf;
                    ahfszz[k] += over_pf*AIZZ[iind][jind][0] * sdens_pf;
                    ahfsxy[k] += over_pf*AIXY[iind][jind][0] * sdens_pf;
                    ahfsxz[k] += over_pf*AIXZ[iind][jind][0] * sdens_pf;
                    ahfsyz[k] += over_pf*AIYZ[iind][jind][0] * sdens_pf;
                  }
		}
              }
            }
          }





	  /*****************************************************
	  *      Computing various properties over a grid      *
	  *****************************************************/
	  
          if (grid)
   	    switch(grid) {

	        /* Electrostatic potential */
	      case 1:
	      for(ix=0;ix<=nix;ix++) {
		pcx = px - (grid_origin[0] + grid_step_x[0]*ix) + grid_step_y[0];
		pcy = py - (grid_origin[1] + grid_step_x[1]*ix) + grid_step_y[1];
		pcz = pz - (grid_origin[2] + grid_step_x[2]*ix) + grid_step_y[2];
	      for(iy=0;iy<=niy;iy++) {
		pcx -= grid_step_y[0];
		pcy -= grid_step_y[1];
                pcz -= grid_step_y[2];
                AI_OSrecurs(pax,pay,paz,pbx,pby,pbz,pcx,pcy,pcz,gamma,iang,jang);
                for(ibf=0;ibf<nbfi;ibf++) {
                  if (eq_shell)
                   jlim = ibf + 1;
                  else
                    jlim = nbfj;
                  for(jbf=0;jbf<jlim;jbf++) {
                    ib = i_1stbf + ibf - 1; jb = j_1stbf + jbf - 1;
		    norm_pf = norm_bf[iang][ibf]*norm_bf[jang][jbf];
		    /*---
		      
		      ATTENTION: Normalization prefactor is included in dens_pf
		      
		     ---*/
                    dens_pf = Ptot[ioff[ib]+jb] * ((ib == jb) ? 1.0 : 2.0) * norm_pf;
                    if (fabs(dens_pf) < PFACCUTOFF) continue;
                    iind = zpow_bf[iang][ibf]*izm + ypow_bf[iang][ibf]*iym +
                           xpow_bf[iang][ibf]*ixm;
                    jind = zpow_bf[jang][jbf]*jzm + ypow_bf[jang][jbf]*jym +
                           xpow_bf[jang][jbf]*jxm;
                    
                    grid_pts[ix][iy] += -over_pf*AI0[iind][jind][0] * dens_pf;
                  }
                }
	      }
	      }
	      break;

		/* Electron density */
	      case 2:
	      if (spin_prop)    /* If spin_prop is set, compute spin density... */
		Density = Pspin;
	      else              /* .. otherwise compute electron density */
		Density = Ptot;
	      for(ix=0;ix<=nix;ix++) {
		pcx = grid_origin[0] + grid_step_x[0]*ix - grid_step_y[0];
		pcy = grid_origin[1] + grid_step_x[1]*ix - grid_step_y[1];
		pcz = grid_origin[2] + grid_step_x[2]*ix - grid_step_y[2];
	      for(iy=0;iy<=niy;iy++) {
		pcx += grid_step_y[0];
		pcy += grid_step_y[1];
		pcz += grid_step_y[2];
		cax = pcx-ax; cbx = pcx-bx;
		cax2 = cax*cax; cbx2 = cbx*cbx;
		cay = pcy-ay; cby = pcy-by;
		cay2 = cay*cay; cby2 = cby*cby;
		caz = pcz-az; cbz = pcz-bz;
		caz2 = caz*caz; cbz2 = cbz*cbz;
		prod_pf = contr[ig]*contr[jg]*exp(-ai*(cax2 + cay2 + caz2))*exp(-bj*(cbx2 + cby2 + cbz2));
		sum1 = 0.0;
		for(ibf=0;ibf<nbfi;ibf++) {
                  if (eq_shell)
                    jlim = ibf + 1;
                  else
                    jlim = nbfj;
                  for(jbf=0;jbf<jlim;jbf++) {
                    ib = i_1stbf + ibf - 1; jb = j_1stbf + jbf - 1;
		    norm_pf = norm_bf[iang][ibf]*norm_bf[jang][jbf];
		    /*---
		      
		      ATTENTION: Normalization prefactor is included in dens_pf
		      
		     ---*/
                    dens_pf = Density[ioff[ib]+jb] * ((ib == jb) ? 1.0 : 2.0) * norm_pf;
		    if (fabs(dens_pf) < PFACCUTOFF) continue;
		    if (eq_shell && ig != jg)
                      dens_pf *= 2.0;
		    /* For **efficiency** reasons I'm trying to avoid calling pow() for up to
		       g-functions through switch() constructs */
                    lx1 = xpow_bf[iang][ibf];
		    switch (lx1) {
		      case 0:
		      cax_l = 1.0;
		      break;
		      case 1:
		      cax_l = cax;
		      break;
		      case 2:
		      cax_l = cax2;
		      break;
		      case 3:
		      cax_l = cax2*cax;
		      break;
		      case 4:
		      cax_l = cax2*cax2;
		      break;
		      default:
		      cax_l = pow(cax,lx1);
		    }
		    ly1 = ypow_bf[iang][ibf];
		    switch (ly1) {
		      case 0:
		      cay_l = 1.0;
		      break;
		      case 1:
		      cay_l = cay;
		      break;
		      case 2:
		      cay_l = cay2;
		      break;
		      case 3:
		      cay_l = cay2*cay;
		      break;
		      case 4:
		      cay_l = cay2*cay2;
		      break;
		      default:
		      cay_l = pow(cay,ly1);
		    }
		    lz1 = zpow_bf[iang][ibf];
		    switch (lz1) {
		      case 0:
		      caz_l = 1.0;
		      break;
		      case 1:
		      caz_l = caz;
		      break;
		      case 2:
		      caz_l = caz2;
		      break;
		      case 3:
		      caz_l = caz2*caz;
		      break;
		      case 4:
		      caz_l = caz2*caz2;
		      break;
		      default:
		      caz_l = pow(caz,lz1);
		    }
		    lx2 = xpow_bf[jang][jbf];
		    switch (lx2) {
		      case 0:
		      cbx_l = 1.0;
		      break;
		      case 1:
		      cbx_l = cbx;
		      break;
		      case 2:
		      cbx_l = cbx2;
		      break;
		      case 3:
		      cbx_l = cbx2*cbx;
		      break;
		      case 4:
		      cbx_l = cbx2*cbx2;
		      break;
		      default:
		      cbx_l = pow(cbx,lx2);
		    }
		    ly2 = ypow_bf[jang][jbf];
		    switch (ly2) {
		      case 0:
		      cby_l = 1.0;
		      break;
		      case 1:
		      cby_l = cby;
		      break;
		      case 2:
		      cby_l = cby2;
		      break;
		      case 3:
		      cby_l = cby2*cby;
		      break;
		      case 4:
		      cby_l = cby2*cby2;
		      break;
		      default:
		      cby_l = pow(cby,ly2);
		    }
		    lz2 = zpow_bf[jang][jbf];
		    switch (lz2) {
		      case 0:
		      cbz_l = 1.0;
		      break;
		      case 1:
		      cbz_l = cbz;
		      break;
		      case 2:
		      cbz_l = cbz2;
		      break;
		      case 3:
		      cbz_l = cbz2*cbz;
		      break;
		      case 4:
		      cbz_l = cbz2*cbz2;
		      break;
		      default:
		      cbz_l = pow(cbz,lz2);
		    }
		    tmp = dens_pf*cax_l*cay_l*caz_l*cbx_l*cby_l*cbz_l;
		    sum1 += tmp;
		  }
		}
		grid_pts[ix][iy] += sum1*prod_pf;
	      }
	      }
	      break;

		/* Electron density gradient vector */
	      case 3:
	      if (spin_prop)    /* If spin_prop is set, compute gradient of spin density... */
		Density = Pspin;
	      else              /* .. otherwise compute gradient of the electron density */
		Density = Ptot;
	      for(ix=0;ix<=nix;ix++) {
	        pcx = grid_origin[0] + grid_step_x[0]*ix - grid_step_y[0];
		pcy = grid_origin[1] + grid_step_x[1]*ix - grid_step_y[1];
		pcz = grid_origin[2] + grid_step_x[2]*ix - grid_step_y[2];
	      for(iy=0;iy<=niy;iy++) {
	        pcx += grid_step_y[0];
		pcy += grid_step_y[1];
		pcz += grid_step_y[2];
		cax = pcx-ax; cbx = pcx-bx;
		cax2 = cax*cax; cbx2 = cbx*cbx;
		cay = pcy-ay; cby = pcy-by;
		cay2 = cay*cay; cby2 = cby*cby;
		caz = pcz-az; cbz = pcz-bz;
		caz2 = caz*caz; cbz2 = cbz*cbz;
		prod_pf = contr[ig]*contr[jg]*exp(-ai*(cax2 + cay2 + caz2))*exp(-bj*(cbx2 + cby2 + cbz2));
		sum1 = sum2 = sum3 = 0.0;
		for(ibf=0;ibf<nbfi;ibf++) {
                  if (eq_shell)
                    jlim = ibf + 1;
                  else
                    jlim = nbfj;
                  for(jbf=0;jbf<jlim;jbf++) {
                    ib = i_1stbf + ibf - 1; jb = j_1stbf + jbf - 1;
		    norm_pf = norm_bf[iang][ibf]*norm_bf[jang][jbf];
		    /*---
		      
		      ATTENTION: Normalization prefactor is included in dens_pf
		      
		     ---*/
                    dens_pf = Density[ioff[ib]+jb] * ((ib == jb) ? 1.0 : 2.0) * norm_pf;
		    if (fabs(dens_pf) < PFACCUTOFF) continue;
		    if (eq_shell && ig != jg)
                      dens_pf *= 2.0;
		    /* For **efficiency** reasons I'm trying to avoid calling pow() for up to
		       g-functions through switch() constructs */
		    lx1 = xpow_bf[iang][ibf];
		    switch (lx1) {
		      case 0:
		      cax_l_1 = cax_l = 1.0;
		      break;
		      case 1:
		      cax_l_1 = 1.0;
		      cax_l = cax;
		      break;
		      case 2:
		      cax_l_1 = cax;
		      cax_l = cax2;
		      break;
		      case 3:
		      cax_l_1 = cax2;
		      cax_l = cax2*cax;
		      break;
		      case 4:
		      cax_l_1 = cax2*cax;
		      cax_l = cax2*cax2;
		      break;
		      default:
		      cax_l_1 = pow(cax,lx1-1);
		      cax_l = cax_l_1*cax;
		    }
		    ly1 = ypow_bf[iang][ibf];
		    switch (ly1) {
		      case 0:
		      cay_l_1 = cay_l = 1.0;
		      break;
		      case 1:
		      cay_l_1 = 1.0;
		      cay_l = cay;
		      break;
		      case 2:
		      cay_l_1 = cay;
		      cay_l = cay2;
		      break;
		      case 3:
		      cay_l_1 = cay2;
		      cay_l = cay2*cay;
		      break;
		      case 4:
		      cay_l_1 = cay2*cay;
		      cay_l = cay2*cay2;
		      break;
		      default:
		      cay_l_1 = pow(cay,ly1-1);
		      cay_l = cay_l_1*cay;
		    }
		    lz1 = zpow_bf[iang][ibf];
		    switch (lz1) {
		      case 0:
		      caz_l_1 = caz_l = 1.0;
		      break;
		      case 1:
		      caz_l_1 = 1.0;
		      caz_l = caz;
		      break;
		      case 2:
		      caz_l_1 = caz;
		      caz_l = caz2;
		      break;
		      case 3:
		      caz_l_1 = caz2;
		      caz_l = caz2*caz;
		      break;
		      case 4:
		      caz_l_1 = caz2*caz;
		      caz_l = caz2*caz2;
		      break;
		      default:
		      caz_l_1 = pow(caz,lz1-1);
		      caz_l = caz_l_1*caz;
		    }
		    lx2 = xpow_bf[jang][jbf];
		    switch (lx2) {
		      case 0:
		      cbx_l_1 = cbx_l = 1.0;
		      break;
		      case 1:
		      cbx_l_1 = 1.0;
		      cbx_l = cbx;
		      break;
		      case 2:
		      cbx_l_1 = cbx;
		      cbx_l = cbx2;
		      break;
		      case 3:
		      cbx_l_1 = cbx2;
		      cbx_l = cbx2*cbx;
		      break;
		      case 4:
		      cbx_l_1 = cbx2*cbx;
		      cbx_l = cbx2*cbx2;
		      break;
		      default:
		      cbx_l_1 = pow(cbx,lx2-1);
		      cbx_l = cbx_l_1*cbx;
		    }
		    ly2 = ypow_bf[jang][jbf];
		    switch (ly2) {
		      case 0:
		      cby_l_1 = cby_l = 1.0;
		      break;
		      case 1:
		      cby_l_1 = 1.0;
		      cby_l = cby;
		      break;
		      case 2:
		      cby_l_1 = cby;
		      cby_l = cby2;
		      break;
		      case 3:
		      cby_l_1 = cby2;
		      cby_l = cby2*cby;
		      break;
		      case 4:
		      cby_l_1 = cby2*cby;
		      cby_l = cby2*cby2;
		      break;
		      default:
		      cby_l_1 = pow(cby,ly2-1);
		      cby_l = cby_l_1*cby;
		    }
		    lz2 = zpow_bf[jang][jbf];
		    switch (lz2) {
		      case 0:
		      cbz_l_1 = cbz_l = 1.0;
		      break;
		      case 1:
		      cbz_l_1 = 1.0;
		      cbz_l = cbz;
		      break;
		      case 2:
		      cbz_l_1 = cbz;
		      cbz_l = cbz2;
		      break;
		      case 3:
		      cbz_l_1 = cbz2;
		      cbz_l = cbz2*cbz;
		      break;
		      case 4:
		      cbz_l_1 = cbz2*cbz;
		      cbz_l = cbz2*cbz2;
		      break;
		      default:
		      cbz_l_1 = pow(cbz,lz2-1);
		      cbz_l = cbz_l_1*cbz;
		    }
		    x2 = cax_l*cbx_l;
		    y2 = cay_l*cby_l;
		    z2 = caz_l*cbz_l;
		      /* X-component of the gradient */
		    sum1 += dens_pf*y2*z2*((lx1*cax_l_1 - 2*ai*cax_l*cax)*cbx_l +
			                   (lx2*cbx_l_1 - 2*bj*cbx_l*cbx)*cax_l );
		      /* Y-component of the gradient */
		    sum2 += dens_pf*x2*z2*((ly1*cay_l_1 - 2*ai*cay_l*cay)*cby_l +
			                   (ly2*cby_l_1 - 2*bj*cby_l*cby)*cay_l );
		      /* Z-component of the gradient */
		    sum3 += dens_pf*x2*y2*((lz1*caz_l_1 - 2*ai*caz_l*caz)*cbz_l +
			                   (lz2*cbz_l_1 - 2*bj*cbz_l*cbz)*caz_l );
		  }
		}
		grid_vecX[ix][iy] += sum1*prod_pf;
		grid_vecY[ix][iy] += sum2*prod_pf;
		grid_vecZ[ix][iy] += sum3*prod_pf;
	      }
	      }
	      break;


	        /* Laplacian of electron density */
	      case 4:
	      if (spin_prop)    /* If spin_prop is set, compute the Laplacian of the spin density... */
		Density = Pspin;
	      else              /* .. otherwise compute the Laplacian of the electron density */
		Density = Ptot;
	      for(ix=0;ix<=nix;ix++) {
	        pcx = grid_origin[0] + grid_step_x[0]*ix - grid_step_y[0];
		pcy = grid_origin[1] + grid_step_x[1]*ix - grid_step_y[1];
		pcz = grid_origin[2] + grid_step_x[2]*ix - grid_step_y[2];
	      for(iy=0;iy<=niy;iy++) {
		pcx += grid_step_y[0];
		pcy += grid_step_y[1];
		pcz += grid_step_y[2];
		cax = pcx-ax; cbx = pcx-bx;
		cax2 = cax*cax; cbx2 = cbx*cbx;
		cay = pcy-ay; cby = pcy-by;
		cay2 = cay*cay; cby2 = cby*cby;
		caz = pcz-az; cbz = pcz-bz;
		caz2 = caz*caz; cbz2 = cbz*cbz;
		prod_pf = contr[ig]*contr[jg]*exp(-ai*(cax2 + cay2 + caz2))*exp(-bj*(cbx2 + cby2 + cbz2));
		sum1 = 0.0;
		for(ibf=0;ibf<nbfi;ibf++) {
                  if (eq_shell)
                    jlim = ibf + 1;
                  else
                    jlim = nbfj;
                  for(jbf=0;jbf<jlim;jbf++) {
                    ib = i_1stbf + ibf - 1; jb = j_1stbf + jbf - 1;
		    norm_pf = norm_bf[iang][ibf]*norm_bf[jang][jbf];
		    /*---
		      
		      ATTENTION: Normalization prefactor is included in dens_pf
		      
		     ---*/
                    dens_pf = Density[ioff[ib]+jb] * ((ib == jb) ? 1.0 : 2.0) * norm_pf;
		    if (fabs(dens_pf) < PFACCUTOFF) continue;
		    if (eq_shell && ig != jg)
                      dens_pf *= 2.0;
		    lx1 = xpow_bf[iang][ibf];
		    ly1 = ypow_bf[iang][ibf];
		    lz1 = zpow_bf[iang][ibf];
		    /* For **efficiency** reasons I'm trying to avoid calling pow() for up to
			 g-functions through switch() constructs */
		    switch (lx1) {
		      case 0:
		      cax_l_2 = cax_l_1 = cax_l = 1.0;
		      break;
		      case 1:
		      cax_l_2 = cax_l_1 = 1.0;
		      cax_l = cax;
		      break;
		      case 2:
		      cax_l_2 = 1.0;
		      cax_l_1 = cax;
		      cax_l = cax2;
		      break;
		      case 3:
		      cax_l_2 = cax;
		      cax_l_1 = cax2;
		      cax_l = cax2*cax;
		      break;
		      case 4:
		      cax_l_2 = cax2;
		      cax_l_1 = cax2*cax;
		      cax_l = cax2*cax2;
		      break;
		      default:
		      cax_l_2 = pow(cax,lx1-2);
		      cax_l_1 = cax_l_2*cax;
		      cax_l = cax_l_1*cax;
		    }
		    switch (ly1) {
		      case 0:
		      cay_l_2 = cay_l_1 = cay_l = 1.0;
		      break;
		      case 1:
		      cay_l_2 = cay_l_1 = 1.0;
		      cay_l = cay;
		      break;
		      case 2:
		      cay_l_2 = 1.0;
		      cay_l_1 = cay;
		      cay_l = cay2;
		      break;
		      case 3:
		      cay_l_2 = cay;
		      cay_l_1 = cay2;
		      cay_l = cay2*cay;
		      break;
		      case 4:
		      cay_l_2 = cay2;
		      cay_l_1 = cay2*cay;
		      cay_l = cay2*cay2;
		      break;
		      default:
		      cay_l_2 = pow(cay,ly1-2);
		      cay_l_1 = cay_l_2*cay;
		      cay_l = cay_l_1*cay;
		    }
		    switch (lz1) {
		      case 0:
		      caz_l_2 = caz_l_1 = caz_l = 1.0;
		      break;
		      case 1:
		      caz_l_2 = caz_l_1 = 1.0;
		      caz_l = caz;
		      break;
		      case 2:
		      caz_l_2 = 1.0;
		      caz_l_1 = caz;
		      caz_l = caz2;
		      break;
		      case 3:
		      caz_l_2 = caz;
		      caz_l_1 = caz2;
		      caz_l = caz2*caz;
		      break;
		      case 4:
		      caz_l_2 = caz2;
		      caz_l_1 = caz2*caz;
		      caz_l = caz2*caz2;
		      break;
		      default:
		      caz_l_2 = pow(caz,lz1-2);
		      caz_l_1 = caz_l_2*caz;
		      caz_l = caz_l_1*caz;
		    }
		    lx2 = xpow_bf[jang][jbf];
		    ly2 = ypow_bf[jang][jbf];
		    lz2 = zpow_bf[jang][jbf];
		    switch (lx2) {
		      case 0:
		      cbx_l_2 = cbx_l_1 = cbx_l = 1.0;
		      break;
		      case 1:
		      cbx_l_2 = cbx_l_1 = 1.0;
		      cbx_l = cbx;
		      break;
		      case 2:
		      cbx_l_2 = 1.0;
		      cbx_l_1 = cbx;
		      cbx_l = cbx2;
		      break;
		      case 3:
		      cbx_l_2 = cbx;
		      cbx_l_1 = cbx2;
		      cbx_l = cbx2*cbx;
		      break;
		      case 4:
		      cbx_l_2 = cbx2;
		      cbx_l_1 = cbx2*cbx;
		      cbx_l = cbx2*cbx2;
		      break;
		      default:
		      cbx_l_2 = pow(cbx,lx2-2);
		      cbx_l_1 = cbx_l_2*cbx;
		      cbx_l = cbx_l_1*cbx;
		    }
		    switch (ly2) {
		      case 0:
		      cby_l_2 = cby_l_1 = cby_l = 1.0;
		      break;
		      case 1:
		      cby_l_2 = cby_l_1 = 1.0;
		      cby_l = cby;
		      break;
		      case 2:
		      cby_l_2 = 1.0;
		      cby_l_1 = cby;
		      cby_l = cby2;
		      break;
		      case 3:
		      cby_l_2 = cby;
		      cby_l_1 = cby2;
		      cby_l = cby2*cby;
		      break;
		      case 4:
		      cby_l_2 = cby2;
		      cby_l_1 = cby2*cby;
		      cby_l = cby2*cby2;
		      break;
		      default:
		      cby_l_2 = pow(cby,ly2-2);
		      cby_l_1 = cby_l_2*cby;
		      cby_l = cby_l_1*cby;
		    }
		    switch (lz2) {
		      case 0:
		      cbz_l_2 = cbz_l_1 = cbz_l = 1.0;
		      break;
		      case 1:
		      cbz_l_2 = cbz_l_1 = 1.0;
		      cbz_l = cbz;
		      break;
		      case 2:
		      cbz_l_2 = 1.0;
		      cbz_l_1 = cbz;
		      cbz_l = cbz2;
		      break;
		      case 3:
		      cbz_l_2 = cbz;
		      cbz_l_1 = cbz2;
		      cbz_l = cbz2*cbz;
		      break;
		      case 4:
		      cbz_l_2 = cbz2;
		      cbz_l_1 = cbz2*cbz;
		      cbz_l = cbz2*cbz2;
		      break;
		      default:
		      cbz_l_2 = pow(cbz,lz2-2);
		      cbz_l_1 = cbz_l_2*cbz;
		      cbz_l = cbz_l_1*cbz;
		    }
		    x2 = cax_l*cbx_l;
		    y2 = cay_l*cby_l;
		    z2 = caz_l*cbz_l;
		      /* d2/dx2 part */
		    sum2 = -y2*z2*
			 ( (lx1*(lx1-1)*cax_l_2 + 4*ai*ai*cax_l*cax2 - 2*ai*(2*lx1+1)*cax_l)*cbx_l +
			   (lx2*(lx2-1)*cbx_l_2 + 4*bj*bj*cbx_l*cbx2 - 2*bj*(2*lx2+1)*cbx_l)*cax_l +
			   2*(lx1*cax_l_1 - 2*ai*cax*cax_l)*(lx2*cbx_l_1 - 2*bj*cbx*cbx_l) );
		      /* d2/dy2 part */
		    sum2 -= x2*z2*
			 ( (ly1*(ly1-1)*cay_l_2 + 4*ai*ai*cay_l*cay2 - 2*ai*(2*ly1+1)*cay_l)*cby_l +
			   (ly2*(ly2-1)*cby_l_2 + 4*bj*bj*cby_l*cby2 - 2*bj*(2*ly2+1)*cby_l)*cay_l +
			   2*(ly1*cay_l_1 - 2*ai*cay*cay_l)*(ly2*cby_l_1 - 2*bj*cby*cby_l) );
		      /* d2/dz2 part */
		    sum2 -= x2*y2*
			 ( (lz1*(lz1-1)*caz_l_2 + 4*ai*ai*caz_l*caz2 - 2*ai*(2*lz1+1)*caz_l)*cbz_l +
			   (lz2*(lz2-1)*cbz_l_2 + 4*bj*bj*cbz_l*cbz2 - 2*bj*(2*lz2+1)*cbz_l)*caz_l +
			   2*(lz1*caz_l_1 - 2*ai*caz*caz_l)*(lz2*cbz_l_1 - 2*bj*cbz*cbz_l) );
		      /* Contract with the density */
		    sum1 += dens_pf*sum2;
		  }
		}
	      grid_pts[ix][iy] += sum1*prod_pf;
	      }
	      }
	      break;
	    }
	  /* End of loops over primitives (ig,jg)*/
        }
      }
      /* End of loops over shells (i,j)*/
    }
  }


  
  if (print_lvl >= PRINTDIPMLEVEL) {
    fprintf(outfile,"\tmu(X) matrix in AO basis :\n");
    print_array(DX_ao,nbfao,outfile);
    fprintf(outfile,"\n\n");
    fprintf(outfile,"\tmu(Y) matrix in AO basis :\n"); 
    print_array(DY_ao,nbfao,outfile);
    fprintf(outfile,"\n\n");
    fprintf(outfile,"\tmu(Z) matrix in AO basis :\n"); 
    print_array(DZ_ao,nbfao,outfile);
    fprintf(outfile,"\n\n");
  }

	/* Computing nuclear contribution to the dipole moment */
  dz_e = dx; dz_e = dx; dz_e = dz;
  for(i=0;i<natom;i++) {
    dx += zvals[i]*geom[i][0];
    dx_n += zvals[i]*geom[i][0];
    dy += zvals[i]*geom[i][1];
    dy_n += zvals[i]*geom[i][1];
    dz += zvals[i]*geom[i][2];
    dz_n += zvals[i]*geom[i][2];
  }
  dtot = sqrt(dx*dx + dy*dy + dz*dz);


	/* Transforming dipole moment integrals to the SO basis */

  if (wrt_dipints) {
    DX_so = init_array(nstri);
    DY_so = init_array(nstri);
    DZ_so = init_array(nstri);
    for(i=0;i<nbfso;i++)
      for(j=0;j<=i;j++) {
        sum1 = 0.0;
        sum2 = 0.0;
        sum3 = 0.0;
        for(k=0;k<nbfao;k++)
          for(l=0;l<nbfao;l++) {
            tmp = usotao[i][k] * usotao[j][l];
            kk = (k > l) ? k : l;
            ll = (k < l) ? k : l;
            sum1 += tmp * DX_ao[ioff[kk]+ll];
            sum2 += tmp * DY_ao[ioff[kk]+ll];
            sum3 += tmp * DZ_ao[ioff[kk]+ll];
          }
        DX_so[ioff[i]+j] = sum1;
        DY_so[ioff[i]+j] = sum2;
        DZ_so[ioff[i]+j] = sum3;
      }
		/* Write them out here for Daniel and dispose */

    rfile(dip_file);
    wwritw(dip_file, (char *) DX_so, sizeof(double)*nstri,0,&junk);
    wwritw(dip_file, (char *) DY_so, sizeof(double)*nstri,junk,&junk);
    wwritw(dip_file, (char *) DZ_so, sizeof(double)*nstri,junk,&junk);
    rclose(dip_file,3);

    if (print_lvl >= PRINTDIPMLEVEL) {
      fprintf(outfile,"\n  mu(X) matrix in SO basis :\n");
      print_array(DX_so,nbfso,outfile);
      fprintf(outfile,"\n  mu(Y) matrix in SO basis :\n");
      print_array(DY_so,nbfso,outfile);
      fprintf(outfile,"\n  mu(Z) matrix in SO basis :\n");
      print_array(DZ_so,nbfso,outfile);
      fprintf(outfile,"\n\n");
    }
    
    free(DX_so);
    free(DY_so);
    free(DZ_so);
    free(DX_ao);
    free(DY_ao);
    free(DZ_ao);
  }

  	/* Computing the total quadrupole moment in traceless tensor form,
  	   and orbital spatial extents */

  if (mpmax > 1) {
    for(i=0;i<natom;i++) {
      x = geom[i][0];
      y = geom[i][1];
      z = geom[i][2];
      mxx += zvals[i] * x * x;
      myy += zvals[i] * y * y;
      mzz += zvals[i] * z * z;
      mxy += zvals[i] * x * y;
      mxz += zvals[i] * x * z;
      myz += zvals[i] * y * z;
    }
    qxx = mxx - (myy+mzz)/2;
    qyy = myy - (mxx+mzz)/2;
    qzz = mzz - (mxx+myy)/2;
    qxy = 1.5 * mxy;
    qxz = 1.5 * mxz;
    qyz = 1.5 * myz;
    
    if (corr) {
      qxxcc = mxxcc - (myycc+mzzcc)/2;
      qyycc = myycc - (mxxcc+mzzcc)/2;
      qzzcc = mzzcc - (mxxcc+myycc)/2;
      qxycc = 1.5 * mxycc;
      qxzcc = 1.5 * mxzcc;
      qyzcc = 1.5 * myzcc;
    }
    
    temp_mat = init_matrix(3,3);
    temp_mat[0][0] = qxx + qxxcc;
    temp_mat[1][1] = qyy + qyycc;
    temp_mat[2][2] = qzz + qzzcc;
    temp_mat[0][1] = temp_mat[1][0] = qxy + qxycc;
    temp_mat[0][2] = temp_mat[2][0] = qxz + qxzcc;
    temp_mat[1][2] = temp_mat[2][1] = qyz + qyzcc;
    sq_rsp(3,3,temp_mat,qvals,1,qvecs,1.0E-14);
    free_matrix(temp_mat,3);
    
	/* Computing orbital spatial extents */
    
    if (read_opdm)	/* If onepdm is read from file - use natural orbitals */
      evec = nmo_ao;
    else
      evec = scf_evec_ao;
    for(i=0;i<nmo;i++) {
      sum1 = sum2 = sum3 = 0.0;
      for(k=0;k<nbfao;k++)
        for(l=0;l<nbfao;l++) {
          tmp = evec[k][i] * evec[l][i];
          kk = (k > l) ? k : l;
          ll = (k < l) ? k : l;
          sum1 += XX[ioff[kk]+ll]*tmp;
          sum2 += YY[ioff[kk]+ll]*tmp;      
          sum3 += ZZ[ioff[kk]+ll]*tmp;      
        }
      MOXX[i] = sum1;
      MOYY[i] = sum2;
      MOZZ[i] = sum3;
    }
    free(XX);
    free(YY);
    free(ZZ);
  }

	/* Computing the total octopole moment in traceless tensor form */

  if (mpmax > 2) {
    for(i=0;i<natom;i++) {
      x=geom[i][0];
      y=geom[i][1];
      z=geom[i][2];
      mxxx += zvals[i]*x*x*x;
      myyy += zvals[i]*y*y*y;
      mzzz += zvals[i]*z*z*z;
      mxxy += zvals[i]*x*x*y;
      mxxz += zvals[i]*x*x*z;
      mxyy += zvals[i]*x*y*y;
      myyz += zvals[i]*y*y*z;
      mxzz += zvals[i]*x*z*z;
      myzz += zvals[i]*y*z*z;
      mxyz += zvals[i]*x*y*z;
    }
    oxxx = mxxx - 1.5*(mxyy+mxzz);
    oyyy = myyy - 1.5*(mxxy+myzz);
    ozzz = mzzz - 1.5*(mxxz+myyz);
    oxxy = 2*mxxy - 0.5*(myyy+myzz);
    oxxz = 2*mxxz - 0.5*(myyz+mzzz);
    oxyy = 2*mxyy - 0.5*(mxxx+mxzz);
    oyyz = 2*myyz - 0.5*(mxxz+mzzz);
    oxzz = 2*mxzz - 0.5*(mxxx+mxyy);
    oyzz = 2*myzz - 0.5*(mxxy+myyy);
    oxyz = 2.5*mxyz;

    if (corr) {
      oxxxcc = mxxxcc - 1.5*(mxyycc+mxzzcc);
      oyyycc = myyycc - 1.5*(mxxycc+myzzcc);
      ozzzcc = mzzzcc - 1.5*(mxxzcc+myyzcc);
      oxxycc = 2*mxxycc - 0.5*(myyycc+myzzcc);
      oxxzcc = 2*mxxzcc - 0.5*(myyzcc+mzzzcc);
      oxyycc = 2*mxyycc - 0.5*(mxxxcc+mxzzcc);
      oyyzcc = 2*myyzcc - 0.5*(mxxzcc+mzzzcc);
      oxzzcc = 2*mxzzcc - 0.5*(mxxxcc+mxyycc);
      oyzzcc = 2*myzzcc - 0.5*(mxxycc+myyycc);
      oxyzcc = 2.5*mxyzcc;
    }
  }


	/* Computing nuclear contribution to the electrostatic properties at the nuclei */

  if (nuc_esp)
    for(k=0;k<natom;k++)
      for(i=0;i<natom;i++)
        if (i != k) {
          x = geom[k][0] - geom[i][0];
          y = geom[k][1] - geom[i][1];
          z = geom[k][2] - geom[i][2];
          r2 = x*x+y*y+z*z;
          r = sqrt(r2);
          phi[k]  += zvals[i]/r;
          ex[k]   += zvals[i]*x/(r*r2);
          ey[k]   += zvals[i]*y/(r*r2);
          ez[k]   += zvals[i]*z/(r*r2);
          dexx[k] += -zvals[i]*(3*x*x-r2)/(r*r2*r2);
          deyy[k] += -zvals[i]*(3*y*y-r2)/(r*r2*r2);
          dezz[k] += -zvals[i]*(3*z*z-r2)/(r*r2*r2);
          dexy[k] += -zvals[i]*3*x*y/(r*r2*r2);
          dexz[k] += -zvals[i]*3*x*z/(r*r2*r2);
          deyz[k] += -zvals[i]*3*x*z/(r*r2*r2);
        }


        /* Computing first-order one-electron relativistic corrections */

    /* Mass-velocity p^4 term */
  massveloc *= -1/(8*_c_au*_c_au);
    /* Darwin term */
  darw = 0.0;
  for(i=0;i<natom;i++)
    darw += zvals[i]*edens[i];
  darw *= M_PI_2/(_c_au*_c_au); 



	/* Computing nuclear contribution to the electrostatic potential over a grid */

  if (grid == 1)
    for(ix=0;ix<=nix;ix++) {
      x = grid_origin[0] + grid_step_x[0]*ix - grid_step_y[0]*iy;
      y = grid_origin[1] + grid_step_x[1]*ix - grid_step_y[1]*iy;
      z = grid_origin[2] + grid_step_x[2]*ix - grid_step_y[2]*iy;
      for(iy=0;iy<=niy;iy++) {
        x += grid_step_y[0]*iy;
        y += grid_step_y[1]*iy;
        z += grid_step_y[2]*iy;
        for(k=0;k<natom;k++) {
          r = sqrt((x-geom[k][0])*(x-geom[k][0]) + 
                   (y-geom[k][1])*(y-geom[k][1]) +
                   (z-geom[k][2])*(z-geom[k][2]));
          if (r > 1.0E-9)
            grid_pts[ix][iy] += zvals[k]/r;
        }
      }
    }

  

        /* Projecting density gradient onto the plane of the grid and computing the magnitude */

  if (grid == 3)
    for(ix=0;ix<=nix;ix++)
      for(iy=0;iy<=niy;iy++) {
	x = grid_vecX[ix][iy]*grid_unit_x[0] + grid_vecY[ix][iy]*grid_unit_x[1] + grid_vecZ[ix][iy]*grid_unit_x[2];
        y = grid_vecX[ix][iy]*grid_unit_y[0] + grid_vecY[ix][iy]*grid_unit_y[1] + grid_vecZ[ix][iy]*grid_unit_y[2];
	grid_pts[ix][iy] = sqrt(x*x + y*y);
	if (edgrad_logscale) {
	  tmp = (double)edgrad_logscale*log(1.0 + grid_pts[ix][iy]);
	  x /= tmp;
	  y /= tmp;
	  grid_pts[ix][iy] /= tmp;
	}
	grid_vecX[ix][iy] = x;
	grid_vecY[ix][iy] = y;
      }


    

	/* Cleaning up */
  
  free_box(MIX,lmax+1,lmax+1);
  free_box(MIY,lmax+1,lmax+1);
  free_box(MIZ,lmax+1,lmax+1);
  if (nuc_esp) {
    free_box(AIX,indmax,indmax);
    free_box(AIY,indmax,indmax);
    free_box(AIZ,indmax,indmax);
    free_box(AIXX,indmax,indmax);
    free_box(AIYY,indmax,indmax);
    free_box(AIZZ,indmax,indmax);
    free_box(AIXY,indmax,indmax);
    free_box(AIXZ,indmax,indmax);
    free_box(AIYZ,indmax,indmax);
  }
  
  if ((grid == 1) || nuc_esp)
    free_box(AI0,indmax,indmax);

}


