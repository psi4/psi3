/*
** PSIFILES.H
**
** This header file contains the definitions of the numbers assigned
**  to various binary files in PSI.  This was created primarily to 
**  help avoid conflicts in the numbering of new PSI files in developmental
**  programs but will grow to encompass some older binary files.
**
** This additional level of abstraction will aid in the maintenance of
**  code.  You are strongly encouraged to refer to files using these
**  definitions rather than the actual numbers; the numbers may change 
**  in the future but the names will not.
**
** Created by C. David Sherrill on 29 April 1998
*/

#define PSIF_SO_INTS        34
#define PSIF_MO_OEI         71
#define PSIF_MO_TEI         72
#define PSIF_MO_OPDM        73
#define PSIF_MO_TPDM        74
#define PSIF_MO_LAG         75
/* PSIF_AO_OPDM also contains AO Lagrangian */
#define PSIF_AO_OPDM        76      
#define PSIF_AO_TPDM        77
#define PSIF_MO_FZC         78

/*--- 
      New additions by Edward Valeev
      Please, contact me if these conflict
      with your definitions
 ---*/
#define PSIF_DSCF           31
#define PSIF_SO_TEI         33
#define PSIF_SO_S           35
#define PSIF_SO_T           36
#define PSIF_SO_V           37
#define PSIF_SO_R12         38
#define PSIF_SO_R12T1       39
#define PSIF_AO_S           40
#define PSIF_AO_MX          41
#define PSIF_AO_MY          42
#define PSIF_AO_MZ          43
#define PSIF_MO_R12         79
#define PSIF_MO_R12T1       80

/*
** Additions for UHF-based transformations.
** -TDC, 6/01
*/
#define PSIF_MO_AA_TEI      81
#define PSIF_MO_BB_TEI      82
#define PSIF_MO_AB_TEI      83
#define PSIF_MO_A_OEI       84
#define PSIF_MO_B_OEI       85
#define PSIF_MO_A_FZC       86
#define PSIF_MO_B_FZC       87

/*
** MO Hessian File (also contains specialized integral and Fock lists.
** See program STABLE for more info.
** -TDC, 7/00
*/
#define PSIF_MO_HESS        88

