#include "scf.h"

#include <libmoinfo/libmoinfo.h>

namespace psi{ namespace mcscf{

void SCF::construct_F()
{
  if(reference == rhf){
    Fc  = H;
    for(int batch = 0; batch < nbatch; ++batch){
      read_Raffanetti("PK",PK,batch);
      construct_G(Dc,G,PK,batch);
      Fc += G;
    }
  }else if(reference == rohf){
    Fc  = H;
    Fo  = H;
    Fo.scale(0.5);
    for(int batch = 0; batch < nbatch; ++batch){
      read_Raffanetti("PK",PK,batch);
      // Dc * PK Contributions
      construct_G(Dc,G,PK,batch);
      Fc += G;
      G.scale(0.5);
      Fo += G;

      // Do * PK Contributions
      construct_G(Do,G,PK,batch,0.5);
      Fc += G;
      G.scale(0.5);
      Fo += G;

      read_Raffanetti("K",K,batch);
      // Do * K Contributions
      construct_G(Do,G,K,batch,0.25);
      Fo += G;
    }
  }else if(reference == tcscf){
    Fc    = H;
    Favg  = H;
    for(int I = 0 ; I < nci; ++I){
      Dsum[I]  = Dc;
      Dsum[I] += Dtc[I];
      Ftc[I] = H;
      Ftc[I].scale(ci[I] * ci[I]);
      H_tcscf[I][I] = 2.0 * dot(Dsum[I],H) + moinfo_scf->get_nuclear_energy();
      for(int J = I + 1; J < nci; ++J)
        H_tcscf[I][J] = H_tcscf[J][I] = 0.0;
    }

    for(int batch = 0; batch < nbatch; ++batch){
      read_Raffanetti("PK",PK,batch);
      // Dc * PK Contributions to the Fock matrices
      construct_G(Dc,G,PK,batch);
      Fc += G;
      for(int I = 0 ; I < nci; ++I){
        T = G;
        T.scale(ci[I] * ci[I]);
        Ftc[I] += T;
      }

      // Dtc * PK Contributions to the Fock matrices
      for(int I = 0 ; I < nci; ++I){
        construct_G(Dtc[I],G,PK,batch,ci[I] * ci[I]);
        Fc += G;
        G.scale(0.5);
        Ftc[I] += G;
      }
       
      // Dsum * PK Contributions to the Hamiltonian
      for(int I = 0 ; I < nci; ++I){
        construct_G(Dsum[I],G,PK,batch);
        H_tcscf[I][I] += dot(Dsum[I],G);

        G.scale(ci[I] * ci[I]);
        Favg += G;
      }

      read_Raffanetti("K",K,batch);
      // Dtc * K Contributions
      for(int I = 0 ; I < nci; ++I){
        construct_G(Dtc[I],G,K,batch);
        T = G;
        T.scale(-0.5 * ci[I] * ci[I]);
        Ftc[I] += T;
        for(int J = 0 ; J < nci; ++J){
          if(I != J){
            T = G;
            T.scale(- ci[I] * ci[J]);
            Ftc[J] += T;
            H_tcscf[I][J] -= dot(Dtc[J],G);
          }
        }
      }

//     // Compute off-diagonal elements of H
//     for(int I = 0 ; I < nci; ++I){
//       for(int J = I + 1; J < nci; ++J){
//         construct_G(Dtc[I],G,"K");
//         H_tcscf[I][J] = H_tcscf[J][I] = - dot(Dtc[J],G);
//       }
//     }

    }
  }
}

// void SCF::construct_F()
// {
//   // Fc
//   construct_G(Dc,G,"PK");
//   Fc  = H;
//   Fc += G;
// 
//   if(reference == rohf){
//     // Dc * PK Contributions
//     Fo  = Fc;
//     Fo.scale(0.5);
// 
//     // Do * PK Contributions
//     construct_G(Do,G,"PK",0.5);
//     Fc += G;
//     G.scale(0.5);
//     Fo += G;
// 
//     // Do * K Contributions
//     construct_G(Do,G,"K",0.25);
//     Fo += G;
//   }
// 
//   if(reference == tcscf){
//     for(int I = 0 ; I < nci; ++I){
//       construct_G(Dtc[I],G,"PK",ci[I]*ci[I]);
//       Fc += G;
//     }
// 
//     // Dc * PK Contributions
//     construct_G(Dc,G,"PK");
//     for(int I = 0 ; I < nci; ++I){
//       Ftc[I]  = H;
//       Ftc[I] += G;
//       Ftc[I].scale(ci[I] * ci[I]);
//     }
// 
//     // Dtc[I] * PK Contributions
//     for(int I = 0 ; I < nci; ++I){
//       construct_G(Dtc[I],G,"PK",0.5 * ci[I] * ci[I]);
//       Ftc[I] += G;
//     }
// 
//     // Dtc[I] * PK Contributions
//     for(int I = 0 ; I < nci; ++I){
//       construct_G(Dtc[I],G,"K",-0.5 * ci[I] * ci[I]);
//       Ftc[I] += G;
//     }
// 
//     for(int I = 0 ; I < nci; ++I){
//       for(int J = 0 ; J < nci; ++J){
//         if(I != J){
//           construct_G(Dtc[J],G,"K", - ci[I] * ci[J]);
//           Ftc[I] += G;
//         }
//       }
//     }
//   }
// }

}} /* End Namespaces */

//   // Fc
//   construct_G(Dc,G,PK);
//   for(int h =0; h < nirreps; ++h){
//     double** g = G->get_block(h);
//     double** f = Fc->get_block(h);
//     double** h0 = H->get_block(h);
//     for(int m = 0; m < sopi[h]; ++m){
//       for(int n = 0; n < sopi[h]; ++n){
//         f[m][n] = h0[m][n] + g[m][n];
//       }
//     }
//   }
//   construct_G(Dtc[0],G,PK);
//   for(int h =0; h < nirreps; ++h){
//     double** g = G->get_block(h);
//     double** f = Fc->get_block(h);
//     for(int m = 0; m < sopi[h]; ++m){
//       for(int n = 0; n < sopi[h]; ++n){
//         f[m][n] += ci[0] * ci[0] * g[m][n];
//       }
//     }
//   }
//   construct_G(Dtc[1],G,PK);
//   for(int h =0; h < nirreps; ++h){
//     double** g = G->get_block(h);
//     double** f = Fc->get_block(h);
//     for(int m = 0; m < sopi[h]; ++m){
//       for(int n = 0; n < sopi[h]; ++n){
//         f[m][n] += ci[1] * ci[1] * g[m][n];
//       }
//     }
//   }
// //   Fc->print();



// void SCF::construct_F_tcscf()
// {

// 
//   // Ftc[0]
//   construct_G(Dc,G,PK);
//   for(int h =0; h < nirreps; ++h){
//     double** g = G->get_block(h);
//     double** f = Ftc[0]->get_block(h);
//     double** h0 = H->get_block(h);
//     for(int m = 0; m < sopi[h]; ++m){
//       for(int n = 0; n < sopi[h]; ++n){
//         f[m][n] = ci[0] * ci[0] * (h0[m][n] + g[m][n]);
//       }
//     }
//   }
//   construct_G(Dtc[0],G,PK);
//   for(int h =0; h < nirreps; ++h){
//     double** g = G->get_block(h);
//     double** f = Ftc[0]->get_block(h);
//     for(int m = 0; m < sopi[h]; ++m){
//       for(int n = 0; n < sopi[h]; ++n){
//         f[m][n] += 0.5 * ci[0] * ci[0] * g[m][n];
//       }
//     }
//   }
//   construct_G(Dtc[0],G,K);
//   for(int h =0; h < nirreps; ++h){
//     double** g = G->get_block(h);
//     double** f = Ftc[0]->get_block(h);
//     for(int m = 0; m < sopi[h]; ++m){
//       for(int n = 0; n < sopi[h]; ++n){
//         f[m][n] += -0.5 * ci[0] * ci[0] * g[m][n];
//       }
//     }
//   }
//   construct_G(Dtc[1],G,K);
//   for(int h =0; h < nirreps; ++h){
//     double** g = G->get_block(h);
//     double** f = Ftc[0]->get_block(h);
//     for(int m = 0; m < sopi[h]; ++m){
//       for(int n = 0; n < sopi[h]; ++n){
//         f[m][n] += - ci[0] * ci[1] * g[m][n];
//       }
//     }
//   }
// 
//   // Ftc[1]
//   construct_G(Dc,G,PK);
//   for(int h =0; h < nirreps; ++h){
//     double** g = G->get_block(h);
//     double** f = Ftc[1]->get_block(h);
//     double** h0 = H->get_block(h);
//     for(int m = 0; m < sopi[h]; ++m){
//       for(int n = 0; n < sopi[h]; ++n){
//         f[m][n] = ci[1] * ci[1] * (h0[m][n] + g[m][n]);
//       }
//     }
//   }
//   construct_G(Dtc[1],G,PK);
//   for(int h =0; h < nirreps; ++h){
//     double** g = G->get_block(h);
//     double** f = Ftc[1]->get_block(h);
//     for(int m = 0; m < sopi[h]; ++m){
//       for(int n = 0; n < sopi[h]; ++n){
//         f[m][n] += 0.5 * ci[1] * ci[1] * g[m][n];
//       }
//     }
//   }
//   construct_G(Dtc[1],G,K);
//   for(int h =0; h < nirreps; ++h){
//     double** g = G->get_block(h);
//     double** f = Ftc[1]->get_block(h);
//     for(int m = 0; m < sopi[h]; ++m){
//       for(int n = 0; n < sopi[h]; ++n){
//         f[m][n] += -0.5 * ci[1] * ci[1] * g[m][n];
//       }
//     }
//   }
//   construct_G(Dtc[0],G,K);
//   for(int h =0; h < nirreps; ++h){
//     double** g = G->get_block(h);
//     double** f = Ftc[1]->get_block(h);
//     for(int m = 0; m < sopi[h]; ++m){
//       for(int n = 0; n < sopi[h]; ++n){
//         f[m][n] += - ci[0] * ci[1] * g[m][n];
//       }
//     }
//   }
// }
