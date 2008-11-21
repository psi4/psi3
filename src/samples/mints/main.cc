#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>

#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libipv1/ip_lib.h>
#include <libiwl/iwl.h>
#include <libqt/qt.h>

#include <libmints/basisset.h>
#include <libmints/onebody.h>
#include <libmints/twobody.h>
#include <libmints/integral.h>
#include <libmints/factory.h>
#include <libmints/symmetry.h>
#include <libmints/wavefunction.h>

using namespace psi;

extern "C" {
    char *gprgid();
}

FILE *infile = NULL,
     *outfile = NULL;
char *psi_file_prefix = NULL;

char *gprgid()
{
    const char *prgid = "MINTS";
    return const_cast<char*>(prgid);
}

std::string to_string(const int val);   // In libmints/matrix.cc

int main (int argc, char * argv[]) 
{
    psi_start(&infile,&outfile,&psi_file_prefix,argc-1,argv+1,0);
    tstart(outfile);
    
    // Required for libmints, allocates and computes the following:
    // ioff, fac, df, bc
    Wavefunction::initialize_singletons();
    
    Ref<PSIO> psio(new PSIO);
    psiopp_ipv1_config(psio.pointer());
    Ref<Chkpt> chkpt(new Chkpt(psio.pointer(), PSIO_OPEN_OLD));
    
    // Create a new matrix factory
    Ref<MatrixFactory> factory(new MatrixFactory);
    
    // Initialize the factory with data from checkpoint file.
    factory->init_with_chkpt(chkpt);
    
    // Initialize the psi3 timer library.
    timer_init();

    // Needed in the examples.
    int nso = chkpt->rd_nso();

    // Some simple examples for programming meeting.
    //  1. Creating and loading overlap integral matrix from a psi file.
    #if 0
    {
        RefMatrix overlap = factory->create_matrix(PSIF_SO_S);

        // A RefMatrix knows how to read itself in from a psi file.
        overlap.load(psio, PSIF_OEI, NULL, nso);

        // Print it out
        overlap.print();

        // Going out of scope releases all memory.
    }
    #endif
    
    //  2. Diagonalizing overlap matrix.
    #if 0
    {
        RefMatrix overlap = factory->create_matrix(PSIF_SO_S);
    
        // A RefMatrix knows how to read itself in from a psi file.
        overlap.load(psio, PSIF_OEI, NULL, nso);
    
        // Print it out
        overlap.print();
    
        // Diagonalize the overlap matrix, first we need storage for eigen-vectors and -values.
        RefMatrix eigenvectors = factory->create_matrix("Eigenvectors");
        RefVector eigenvalues  = factory->create_vector();
    
        overlap.diagonalize(eigenvectors, eigenvalues);
    
        // Print out the information.
        eigenvectors.eivprint(eigenvalues);

        // Going out of scope releases all memory.
    }
    #endif
    
    //   3. Adding two matrices together
    #if 0
    {
        // Load in kinetic and potential integrals to form one-electron Hamiltonian.
        RefMatrix kinetic = factory->create_matrix(PSIF_SO_T);
        RefMatrix potential = factory->create_matrix(PSIF_SO_V);
    
        kinetic.load(psio, PSIF_OEI, NULL, nso);
        potential.load(psio, PSIF_OEI, NULL, nso);
    
        // Form one-electron Hamiltonian
        RefMatrix H = factory->create_matrix("One-electron Hamiltonian");
        H.copy(kinetic + potential);
    
        // Print out the matrix
        H.print();
    
        // Going out of scope releases all memory.
    }
    #endif
    
    // This code block is for testing the basis set and one electron integral codes.
    //#if 0
    {
        // Create a basis set object and have it initialize itself using the checkpoint file
        Ref<BasisSet> basis(new BasisSet(chkpt));
            
        char *label = chkpt->rd_label();
        fprintf(outfile, "\nLabel: %s\n\n", label);
        free(label);
    
        // Print the molecule
        basis->molecule()->print();
        int natom = basis->molecule()->natom();
    
        // Initialize an integral object
        Ref<IntegralFactory> integral(new IntegralFactory(basis, basis, basis, basis));
        Ref<OneBodyInt> deriv = integral->dipole(1);
        
        RefSimpleMatrixArray in = new RefSimpleMatrix[3];
        RefSimpleMatrixArray o1 = new RefSimpleMatrix[3*3*natom];
        
        for (int i=0; i<3; ++i) {
            std::string name = "SO-basis " + to_string(i);
            in[i] = factory->create_simple_matrix(name);
        }
        for (int i=0; i<3*3*natom; ++i) {
            std::string name = "SO-basis Derivative " + to_string(i);
            o1[i] = factory->create_simple_matrix(name);
        }
    
        deriv->compute(in);
        deriv->compute_deriv1(o1);
        
        for (int i=0; i<3; ++i)
            in[i]->print();
        fprintf(outfile, "Derivatives:\n");
        for (int i=0; i<3*3*natom; ++i)
            o1[i]->print();        
    }
    //#endif
    
    // This code block is for running MP2
    #if 0
    {
        // Create a basis set object and have it initialize itself using the checkpoint file
        Ref<BasisSet> basis(new BasisSet(chkpt));
        // Ref<Symmetry> symmetry(new Symmetry(chkpt));
        
        char *label = chkpt->rd_label();
        fprintf(outfile, "Label from file32: %s\n\n", label);
        free(label);
    
        // Needed in the example.
        double escf = chkpt->rd_escf();
        
        int nirrep = chkpt->rd_nirreps();
        if (nirrep != 1) {
            fprintf(stderr, "Please run in C1 symmetry.\n");
            return EXIT_FAILURE;
        }
        
        int nocc = 0;
        ip_data(const_cast<char*>("DOCC"), const_cast<char*>("%d"), &nocc, 1, 0);
        if (nocc == 0) {
            fprintf(stderr, "Missing DOCC in input file.\n");
            return EXIT_FAILURE;
        }
        
        // Print the molecule
        basis->molecule()->print();
    
        // Initialize an integral object
        Ref<IntegralFactory> integral(new IntegralFactory(basis, basis, basis, basis));
        Ref<TwoBodyInt> eri = integral->eri();

        RefMatrix C = factory->create_matrix("MO coefficients");
                
        // Compute MP2 in C1 symmetry
        double **vectors = chkpt->rd_scf();
        if (vectors == NULL) {
            fprintf(stderr, "Could not find MO coefficients. Run cscf first.\n");
            return EXIT_FAILURE;
        }
        C.set(const_cast<const double**>(vectors));
        free_block(vectors);
        
        // Storage for the AOs
        double *pqrs = new double[nao*nao*nao*nao];
        memset(pqrs, 0, sizeof(double)*nao*nao*nao*nao);
        const double *buffer = eri->buffer();
        
        fprintf(outfile, "  Computing integrals..."); fflush(outfile);
        
        int nshell = basis->nshell();
        for (int P = 0; P < nshell; P++) {
            int nump = basis->shell(P)->nfunction();
            
            for (int Q = 0; Q < nshell; ++Q) {
                int numq = basis->shell(Q)->nfunction();
                
                for (int R = 0; R < nshell; ++R) {
                    int numr = basis->shell(R)->nfunction();
                    
                    for (int S = 0; S < nshell; ++S) {
                        int nums = basis->shell(S)->nfunction();
                        
                        eri->compute_shell(P, Q, R, S);
                        
                        int index = 0;
                        for(int p=0; p < nump; p++) {
                            int op = basis->shell(P)->function_index()+p;
                            
                            for(int q = 0; q < numq; q++) {
                                int oq = basis->shell(Q)->function_index()+q;
                                
                                for(int r = 0; r < numr; r++) {
                                    int oor = basis->shell(R)->function_index()+r;
                                    
                                    for(int s = 0; s < nums; s++,index++) {
                                        int os = basis->shell(S)->function_index()+s;
                                        
                                        // printf("op = %d oq = %d oor =%d os = %d\n", op, oq, oor, os);
                                        
                                        int ipqrs = (((op*nao+oq)*nao+oor)*nao+os);
                                        pqrs[ipqrs] = buffer[index];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        fclose(ints_out);
        fprintf(outfile, "done.\nUsed %d memory.\n", nao*nao*nao*nao); fflush(outfile);
        
        // Needed in MP2.
         int nmo = chkpt->rd_nmo();
         
         // Transform from the AO to MO basis.
         double *ijkl = new double[nmo*nmo*nmo*nmo];
         int idx = 0;
         fprintf(outfile, "  Transforming (N^8) integrals..."); fflush(outfile);
         for(int i = 0; i < nmo; i++) {
             for(int j = 0; j < nmo; j++) {
                 for(int k = 0; k < nmo; k++) {
                     for(int l = 0; l < nmo; l++, idx++) {
                         
                         ijkl[idx] = 0.0;
                         
                         int index = 0;
                         for(int p = 0; p < nao; p++) {
                             for(int q = 0; q < nao; q++) {
                                 for(int r = 0; r < nao; r++) {
                                     for(int s = 0; s < nao; s++,index++) {
                                         
                                         ijkl[idx] += C->get(0, p, i) * C->get(0, q, j) * C->get(0, r, k) * C->get(0, s, l) * pqrs[index];
                                     }
                                 }
                             }
                         }
                         
                     }
                 }
             }
             fprintf(outfile, "%d%%...", (int)((double)i / (double)nmo * 100.0)); fflush(outfile);
         }
         fprintf(outfile, "done.\n"); fflush(outfile);
        
         // Done with AO integrals
         delete[] pqrs;
        
         // Get the evals:
         double *evals = chkpt->rd_evals();
                 
         fprintf(outfile, "  Computing MP2 energy..."); fflush(outfile);
         double energy = 0.0;
         for(int i=0; i < nocc; i++) {
             for(int j=0; j < nocc; j++) {
                 for(int a=nocc; a < nmo; a++) {
                     for(int b=nocc; b < nmo; b++) {
                         
                         int iajb = (((i*nmo+a)*nmo+j)*nmo+b);
                         int ibja = (((i*nmo+b)*nmo+j)*nmo+a);
                         
                         energy += (2.0 * ijkl[iajb] - ijkl[ibja]) * ijkl[iajb]/
                         (evals[i] + evals[j] - evals[a] - evals[b]);
                         
                     }
                 }
             }
         }
         fprintf(outfile, "done.\n\n"); fflush(outfile);
         
         // Done with everything.
         delete[] evals;
         delete[] ijkl;
         
         fprintf(outfile, "  SCF energy:             %20.14f\n", escf);
         fprintf(outfile, "  MP2 correlation energy: %20.14f\n", energy);
         fprintf(outfile, "  Total MP2 energy:       %20.14f\n", escf + energy);
    }
    #endif
    
    // This block of code is testing ERI derivatives
    #if 0
    {
        // Create a basis set object and have it initialize itself using the checkpoint file
        Ref<BasisSet> basis(new BasisSet(chkpt));

        char *label = chkpt->rd_label();
        fprintf(outfile, "Label from file32: %s\n\n", label);
        free(label);

        // Needed in the example.
        int nao = chkpt->rd_nso();

        int nirrep = chkpt->rd_nirreps();
        if (nirrep != 1) {
            fprintf(stderr, "Please run in C1 symmetry.\n");
            return EXIT_FAILURE;
        }

        // Print the molecule
        basis->molecule()->print();
        int natom = basis->molecule()->natom();
        //basis->print();

        // Initialize an integral object
        Ref<IntegralFactory> integral(new IntegralFactory(basis, basis, basis, basis));
        Ref<TwoBodyInt> eri = integral->eri_deriv();

        const double *buffer = eri->buffer();

        fprintf(outfile, "  Computing integral derivatives..."); fflush(outfile);

        int nshell = basis->nshell();
        FILE *ints_out = fopen("matrix.integral_derivatives", "w");

        for (int P = 0; P < nshell; P++) {
            int nump = basis->shell(P)->nfunction();

            for (int Q = 0; Q < nshell; ++Q) {
                int numq = basis->shell(Q)->nfunction();

                for (int R = 0; R < nshell; ++R) {
                    int numr = basis->shell(R)->nfunction();

                    for (int S = 0; S < nshell; ++S) {
                        int nums = basis->shell(S)->nfunction();

                        eri->compute_shell(P, Q, R, S);

                        size_t size = nump * numq * numr * nums;
                        
                        int index = 0;
                        for(int p=0; p < nump; p++) {
                            int op = basis->shell(P)->function_index()+p;
                        
                            for(int q = 0; q < numq; q++) {
                                int oq = basis->shell(Q)->function_index()+q;
                        
                                for(int r = 0; r < numr; r++) {
                                    int oor = basis->shell(R)->function_index()+r;
                        
                                    for(int s = 0; s < nums; s++,index++) {
                                        int os = basis->shell(S)->function_index()+s;
                        
                                        for (int i = 0; i < 3*natom; ++i)
                                            if (fabs(buffer[i*size+index]) > 1.0e-14)
                                                fprintf(ints_out, " A%d  %3d %3d %3d %3d %20.14f\n", i, op, oq, oor, os, buffer[i*size+index]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        fclose(ints_out);
        fprintf(outfile, "done.\n"); fflush(outfile);
    }
    #endif
    
    // Shut down psi. 
    timer_done();
    tstop(outfile);
    psi_stop(infile, outfile, psi_file_prefix);
    
    return EXIT_SUCCESS;
}
