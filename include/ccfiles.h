/* The miscellaneous CC information file */
#define CC_INFO        100
/* One-electron integral files */
#define CC_OEI         101
/* Two-electron integral files */   /* pqnum  rsnum */
#define CC_AINTS       102
#define CC_BINTS       103
#define CC_CINTS       104
#define CC_DINTS       105
#define CC_EINTS       106
#define CC_FINTS       107

/* Two-electron amplitudes, intermediates, and densities */
#define CC_DENOM       108
#define CC_TAMPS       109
#define CC_GAMMA       110
#define CC_MISC        111
#define CC_HBAR        112

#define CC_OEI_NEW     113
#define CC_GAMMA_NEW   114
#define CC_AINTS_NEW   115
#define CC_BINTS_NEW   116
#define CC_CINTS_NEW   117
#define CC_DINTS_NEW   118
#define CC_EINTS_NEW   119
#define CC_FINTS_NEW   120

/* ground state lambda and intermediates for excited states */
#define CC_LAMBDA      121

/* converged eigenvectors of hbar */
#define CC_RAMPS       122
#define CC_LAMPS       123

#define CC_LR          124

#define CC_DIIS_ERR    125
#define CC_DIIS_AMP    126

#define CC_TMP         127
#define CC_TMP0        128
#define CC_TMP1        129
#define CC_TMP2        130
#define CC_TMP3        131
#define CC_TMP4        132
#define CC_TMP5        133
#define CC_TMP6        134
#define CC_TMP7        135
#define CC_TMP8        135
#define CC_TMP9        137
#define CC_TMP10       138
#define CC_TMP11       139
/* temporary files for CCEOM and CCLAMBDA */
#define EOM_D          140
#define EOM_CME        141
#define EOM_Cme        142
#define EOM_CMNEF      143
#define EOM_Cmnef      144
#define EOM_CMnEf      145
#define EOM_SIA        146
#define EOM_Sia        147
#define EOM_SIJAB      148
#define EOM_Sijab      149
#define EOM_SIjAb      150
#define EOM_TMP        151
#define EOM_R          152

/* Markers for the first and last file numbers */
#define CC_MIN  CC_INFO
#define CC_MAX  EOM_R
