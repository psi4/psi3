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
#define CC_LAMPS       110
#define CC_HBAR        111
#define CC_GAMMA       112
#define CC_MISC        113

/* Some temporary files for short-lived intermediates --- I'll eventually
   pare this down to only one or two */
#define CC_TMP         114
#define CC_TMP0        115
#define CC_TMP1        116
#define CC_TMP2        117
#define CC_TMP3        118
#define CC_TMP4        119
#define CC_TMP5        120
#define CC_TMP6        121
#define CC_TMP7        122
#define CC_TMP8        123
#define CC_TMP9        124
#define CC_TMP10       125
#define CC_TMP11       126
#define CC_OEI_NEW     127
#define CC_GAMMA_NEW   128
#define CC_AINTS_NEW   129
#define CC_BINTS_NEW   130
#define CC_CINTS_NEW   131
#define CC_DINTS_NEW   132
#define CC_EINTS_NEW   133
#define CC_FINTS_NEW   134

/* Markers for the first and last file numbers */
#define CC_MIN  CC_INFO
#define CC_MAX  CC_FINTS_NEW
