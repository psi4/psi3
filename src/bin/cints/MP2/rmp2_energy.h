#define LOCK_RS_SHELL 0      /*--- When updating (js|ia) and (jr|ia) in (JS|IA) lock blocks corresponding to the entire shell
			       blocks or just the appropriate basis functions (more fine-grained in the second case) ---*/

void rmp2_energy(void);

typedef struct {
    int num_i_per_ibatch;
    int num_ibatch;
    int num_arrived;
    double Emp2_0, Emp2_1;
    double *emp2_0, *emp2_1;
} RMP2_Status_t;
