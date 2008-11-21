#include <cstdlib>
#include <cmath>
#include <libmints/gshell.h>

#include <libmints/wavefunction.h>

#include <psiconfig.h>

using namespace psi;

extern FILE *outfile;

double norm(int x1, int x2,double c,double ss)
{
    if (x1 < x2) 
        return norm(x2,x1,c,ss);
    if (x1 == 1) {
        if (x2 == 1)
            return c * ss;
        else
            return 0.0;
    }
    if (x1 == 0)
        return ss;
    return c * ( (x1-1) * norm(x1-2,x2,c,ss) + (x2 * norm(x1-1,x2-1,c,ss)));
}


GaussianShell::GaussianShell(int ncn, int nprm, double* e, int* am, GaussianType pure,
    double** c, int nc, Vector3& center, int start, PrimitiveType pt):
    nprimitives_(nprm), ncontractions_(ncn), nc_(nc), center_(center), start_(start), sym_transfrom_(0)
{
    puream_ = new int[ncontraction()];
    for (int i=0; i<ncontraction(); ++i) {
        puream_[i] = (pure == Pure);
    }
 
    copy_data(am, e, c);
    // Compute the number of basis functions in this shell
    init_data();
    
    // Convert the coefficients to coefficients for unnormalized primitives, if needed
    // if (pt == Normalized)
    //     convert_coefficients();
        
    // Compute the normalization constants
    if (pt == Unnormalized)
        normalize_shell();
}

GaussianShell::~GaussianShell()
{
    delete[] l_;
    delete[] puream_;
    delete[] exp_;
    
    for (int i=0; i<ncontractions_; ++i)
        delete[] coef_[i];
    
    delete[] coef_;
    
    if (sym_transfrom_)
        delete[] sym_transfrom_;
}

// expects coef to be in primitive x contraction format. the data is transposed here
void GaussianShell::copy_data(int *l, double *exp, double **coef)
{
    l_ = new int[ncontraction()];
    for (int c=0; c<ncontraction(); ++c)
        l_[c] = l[c];
    
    exp_ = new double[nprimitive()];
    coef_ = new double*[ncontraction()];
    for (int p=0; p<nprimitive(); ++p) {
        exp_[p] = exp[p];
    }
    for (int c=0; c<ncontraction(); ++c) {
        coef_[c] = new double[nprimitive()];
        for (int p=0; p<nprimitive(); ++p) {
            // Yes, I want it stored c x p, but it came in from chkpt as p x c
            coef_[c][p] = coef[p][c];
        }
    }
}

void GaussianShell::convert_coefficients()
{
    int i,gc;
    double c,ss;

    // Convert the contraction coefficients from coefficients over
    // normalized primitives to coefficients over unnormalized primitives
    for (gc=0; gc<ncontractions_; gc++) {
        for (i=0; i<nprimitives_; i++) {
            c = 0.25/exp_[i];
            ss = pow(M_PI/(exp_[i]+exp_[i]),1.5);
            coef_[gc][i] *= 1.0/sqrt(::norm(l_[gc],l_[gc],c,ss));
        }
    }
}

double GaussianShell::shell_normalization(int gs)
{
    int i,j;
    double result,c,ss;

    result = 0.0;
    for (i=0; i<nprimitives_; i++) {
        for (j=0; j<nprimitives_; j++) {
            c = 0.50/(exp_[i] + exp_[j]);
            ss = pow(M_PI/(exp_[i]+exp_[j]),1.5);
            result += coef_[gs][i] * coef_[gs][j] *
                ::norm(l_[gs],l_[gs],c,ss);
        }
    }

    return 1.0/sqrt(result);
}

void GaussianShell::normalize_shell()
{
    int i, gs;
    
    for (gs = 0; gs < ncontractions_; ++gs) {
        double normalization = shell_normalization(gs);
        for (i = 0; i < nprimitives_; ++i) {
            coef_[gs][i] *= normalization;
            #ifdef DEBUG
            fprintf(outfile, "New c = %20.16f for %20.16f\n", coef_[gs][i], exp_[i]);
            #endif
        }
    }
}

int GaussianShell::nfunction(int c) const
{
    return INT_NFUNC(puream_[c], l_[c]);
}

void GaussianShell::init_data()
{
    int max = 0;
    int min = 0;
    int nc = 0;
    int nf = 0;
    has_pure_ = false;
    
    for (int i=0; i<ncontraction(); ++i) {
        int maxi = l_[i];
        if (max < maxi)
            max = maxi;
            
        int mini = l_[i];
        if (min > mini || i == 0)
            min = mini;
            
        nc += ncartesian(i);
        nf += nfunction(i);
        
        if (is_pure(i))
            has_pure_ = true;
    }
    
    max_am_ = max;
    min_am_ = min;
    ncartesians_ = nc;
    nfunctions_ = nf;
}

void GaussianShell::print(FILE *out) const
{
    fprintf(out, "      Number of contractions: %d\n", ncontraction());
    fprintf(out, "      Number of primitives: %d\n", nprimitive());
    fprintf(out, "      Number of Cartesian Gaussians: %d\n", ncartesian());
    fprintf(out, "      Spherical Harmonics?: %s\n", has_pure() ? "true" : "false");
    if (max_am() == min_am())
        fprintf(out, "      Angular momentum: %d\n", max_am());
    else {
        fprintf(out, "      Max angular momentum: %d\n", max_am());
        fprintf(out, "      Min angular momentum: %d\n", min_am());
    }
    fprintf(out, "      Exponent       ");
    for(int c=0; c<ncontraction(); c++)
        fprintf(out, " Contr. %3d  ",c);
    fprintf(out, "\n");
    for(int p=0; p<nprimitive(); p++) {
        fprintf(out, "      %15.10f",exp_[p]);
        for(int c=0; c<ncontraction(); c++)
            fprintf(out, " %12.9f",coef_[c][p]);
        fprintf(out, "\n");
    }
    fprintf(out, "\n");
}

double GaussianShell::normalize(int l, int m, int n)
{
    static int use_cca_integrals_standard = (PSI_INTEGRALS_STANDARD == 1);
    if (use_cca_integrals_standard) {
        return 1.0;
    } else {
        double numer = df[2*l_[0]];
        double denom = df[2*l] * df[2*m] * df[2*n];
        // printf("l_=%d, l=%d, m=%d, n=%d, norm=%15.10f\n", l_[0], l, m, n, sqrt(numer/denom));
        return sqrt(numer/denom);        
    }
}

void GaussianShell::set_sym_transform(int nirreps, int *vec)
{
    if (sym_transfrom_)
        delete[] sym_transfrom_;
    
    sym_transfrom_ = new int[nirreps];
    for (int i=0; i<nirreps; ++i)
        sym_transfrom_[i] = vec[0];
}
