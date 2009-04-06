#ifndef _psi_src_lib_libmints_matrix_h_
#define _psi_src_lib_libmints_matrix_h_

/*!
    \file libmints/matrix.h
    \ingroup MINTS
*/

#include <cstdio>
#include <string>
#include <cstring>

#include <libmints/vector.h>
#include <libmints/ref.h>

#include <libpsio/psio.hpp>

extern FILE *outfile;

namespace psi {

class MatrixFactory;
class SimpleMatrix;
class RefSimpleMatrix;

//! Makes using matrices just a little earlier.
//! Using a matrix factory makes creating these a breeze.
class Matrix {
protected:
    /// Matrix data
    double ***matrix_;
    /// Number of irreps
    int nirreps_;
    /// Rows per irrep array
    int *rowspi_;
    /// Columns per irrep array
    int *colspi_;
    /// Name of the matrix
    std::string name_;

    /// Allocates matrix_
    void alloc();
    /// Release matrix_
    void release();

    /// Copies data from the passed matrix to this matrix_
    void copy_from(double ***);

    /// allocate a block matrix -- analogous to libciomr's block_matrix
    static double** matrix(int nrow, int ncol) {
        double** mat = (double**) malloc(sizeof(double*)*nrow);
        const size_t size = sizeof(double)*nrow*ncol;
        mat[0] = (double*) malloc(size);
        //bzero((void*)mat[0],size);
        memset((void *)mat[0], '\0', size);
        for(int r=1; r<nrow; ++r) mat[r] = mat[r-1] + ncol;
        return mat;
    }
    /// free a (block) matrix -- analogous to libciomr's free_block
    static void free(double** Block) {
        ::free(Block[0]);  ::free(Block);
    }

    void print_mat(double **a, int m, int n, FILE *out);

public:
    /// Default constructor, zeros everything out
    Matrix();
    /// Constructor, zeros everything out, sets name_
    Matrix(std::string name);
    /// Explicit copy reference constructor
    explicit Matrix(const Matrix& copy);
    /// Explicit copy pointer constructor
    explicit Matrix(const Matrix* copy);
    /// Constructor, sets up the matrix
    Matrix(int nirreps, int *rowspi, int *colspi);
    /// Constructor, sets name_, and sets up the matrix
    Matrix(std::string name, int nirreps, int *rowspi, int *colspi);

    /// Destructor, frees memory
    ~Matrix();

    /// Creates an exact copy of the matrix and returns it.
    Matrix* clone() const;
    /// Copies cp's data onto this
    void copy(Matrix* cp);

    /// Load a matrix from a PSIO object from fileno with tocentry of size nso
    bool load(Ref<psi::PSIO>& psio, unsigned int fileno, char *tocentry, int nso);

    /// Saves the matrix in ASCII format to filename
    void save(const char *filename, bool append=true, bool saveLowerTriangle = true, bool saveSubBlocks=false);
    void save(std::string filename, bool append=true, bool saveLowerTriangle = true, bool saveSubBlocks=false) {
        save(filename.c_str(), append, saveSubBlocks);
    }
    /// Saves the block matrix to PSIO object with fileno and with the toc position of the name of the matrix
    void save(Ref<psi::PSIO>& psio, unsigned int fileno, bool saveSubBlocks=true);

    /// Set every element of matrix_ to val
    void set(double val);
    /// Copies lower triangle tri to matrix_, calls tri_to_sq
    void set(const double *tri);
    /// Copies sq to matrix_
    void set(const double **sq);
    /// Set a single element of matrix_
    void set(int h, int m, int n, double val) { matrix_[h][m][n] = val; }
    /// Set the diagonal of matrix_ to vec
    void set(Vector *vec);
    /// Returns a single element of matrix_
    double get(int h, int m, int n) { return matrix_[h][m][n]; }
    /// Returns matrix_
    double **to_block_matrix() const;
    /// Converts this to a full non-symmetry-block matrix
    SimpleMatrix *to_simple_matrix();

    /// Sets the name of the matrix, used in print(...) and save(...)
    void set_name(std::string name) {
        name_ = name;
    };

    /// Print the matrix using print_mat
    void print(FILE *out = outfile, char *extra=NULL);
    /// Print the matrix with corresponding eigenvalues below each column
    void eivprint(Vector *values, FILE *out = outfile);
    /// Returns the rows per irrep array
    int *rowspi() const {
        return rowspi_;
    }
    /// Returns the columns per irrep array
    int *colspi() const {
        return colspi_;
    }
    /// Returns the number of irreps
    int nirreps() const {
        return nirreps_;
    }
    /// Set this to identity
    void set_to_identity();
    /// Zeros this out
    void zero();
    /// Zeros the diagonal
    void zero_diagonal();

    // Math routines
    /// Returns the trace of this
    double trace();
    /// Creates a new matrix which is the transpose of this
    Matrix *transpose();

    /// Adds a matrix to this
    void add(const Matrix*);
    /// Subtracts a matrix from this
    void subtract(const Matrix*);
    /// Multiplies the two arguments and adds their result to this
    void accumulate_product(const Matrix*, const Matrix*);
    /// Scales this matrix
    void scale(double);
    /// Returns the sum of the squares of this
    double sum_of_squares();
    /// Add val to an element of this
    void add(int h, int m, int n, double val) { matrix_[h][m][n] += val; }
    /// Scale row m of irrep h by a
    void scale_row(int h, int m, double a);
    /// Scale column n of irrep h by a
    void scale_column(int h, int n, double a);
    /// Transform a by transformer save result to this
    void transform(Matrix* a, Matrix* transformer);
    /// Transform this by transformer
    void transform(Matrix* transformer);
    /// Back transform a by transformer save result to this
    void back_transform(Matrix* a, Matrix* transformer);
    /// Back transform this by transformer
    void back_transform(Matrix* transformer);

    /// Returns the vector dot product of this by rhs
    double vector_dot(Matrix* rhs);

    /// General matrix multiply, saves result to this
    void gemm(bool transa, bool transb, double alpha, const Matrix* a, const Matrix* b, double beta);
    /// Diagonalize this places eigvectors and eigvalues must be created by caller.
    void diagonalize(Matrix* eigvectors, Vector* eigvalues);
};

//! Matrix reference wrapped class.
class RefMatrix : public Ref< Matrix > {
public:
    RefMatrix();
    RefMatrix(Matrix *o);
    RefMatrix(const RefMatrix& o);

    Matrix* clone() const;
    void copy(const RefMatrix& cp);

    bool load(Ref<psi::PSIO>& psio, unsigned int fileno, char *tocentry, int nso);

    void save(const char *filename, bool append=true, bool saveLowerTriangle = true, bool saveSubBlocks=false);
    void save(std::string filename, bool append=true, bool saveLowerTriangle = true, bool saveSubBlocks=false);
    void save(Ref<psi::PSIO>& psio, unsigned int fileno, bool saveSubBlocks=true);

    void set(double val);
    void set(const double *tri);
    void set(const double **sq);
    void set(int h, int m, int n, double val);
    void set(const RefVector& vec);
    double get(int h, int m, int n);
    double **to_block_matrix();
    RefSimpleMatrix to_simple_matrix();

    void set_name(std::string name);

    void print(FILE *out = outfile, char *extra=NULL);
    void eivprint(RefVector& values, FILE *out=outfile);
    int *rowspi() const;
    int *colspi() const;
    int nirreps() const;
    void set_to_identity();
    void zero();
    void zero_diagonal();

    // Math routines
    double trace();
    RefMatrix transpose();

    void add(const RefMatrix&);
    void subtract(const RefMatrix&);
    void scale(double);
    double sum_of_squares();
    void add(int h, int m, int n, double val);
    void scale_row(int h, int m, double a);
    void scale_column(int h, int n, double a);
    void transform(RefMatrix& a, RefMatrix& transformer);
    void transform(RefMatrix& transformer);
    void back_transform(RefMatrix& a, RefMatrix& transformer);
    void back_transform(RefMatrix& transformer);

    double vector_dot(RefMatrix& rhs);

    void gemm(bool transa, bool transb, double alpha, const RefMatrix& a, const RefMatrix& b, double beta);
    void diagonalize(RefMatrix& eigvectors, RefVector& eigvalues);

    // Make this refer to m
    RefMatrix& operator=(Matrix* m);
    // Make this and m refer to the same matrix
    RefMatrix& operator=(const RefMatrix& m);

    // Multiply this by a matrix and return a matrix
    RefMatrix operator*(const RefMatrix&) const;
    // Multiply this by a scalar and return the result
    RefMatrix operator*(double) const;
    // Matrix addition
    RefMatrix operator+(const RefMatrix&) const;
    // Matrix subtraction
    RefMatrix operator-(const RefMatrix&) const;
};

//! Simple matrix class. Not symmetry blocked.
class SimpleMatrix
{
protected:
    /// Matrix data
    double **matrix_;
    /// Number of rows and columns
    int rows_, cols_;
    /// Nae of the matrix
    std::string name_;

    /// Allocates matrix_
    void alloc();
    /// Releases matrix_
    void release();

    /// Copies data from the passed matrix to this matrix_
    void copy_from(double **);

    /// allocate a block matrix -- analogous to libciomr's block_matrix
    static double** matrix(int nrow, int ncol) {
        double** mat = (double**) malloc(sizeof(double*)*nrow);
        const size_t size = sizeof(double)*nrow*ncol;
        mat[0] = (double*) malloc(size);
        //bzero((void*)mat[0],size);
        memset((void *)mat[0], '\0', size);
        for(int r=1; r<nrow; ++r) mat[r] = mat[r-1] + ncol;
        return mat;
    }
    /// free a (block) matrix -- analogous to libciomr's free_block
    static void free(double** Block) {
        ::free(Block[0]);  ::free(Block);
    }

public:
    /// Default constructor, zeros everything out
    SimpleMatrix();
    /// Constructor, zeros everything out, sets name_
    SimpleMatrix(std::string name);
    /// Explicit copy reference constructor
    explicit SimpleMatrix(const SimpleMatrix& copy);
    /// Explicit copy pointer constructor
    explicit SimpleMatrix(const SimpleMatrix* copy);
    /// Constructor, sets up the matrix
    SimpleMatrix(int rows, int cols);
    /// Constructor, sets name_, and sets up the matrix
    SimpleMatrix(std::string name, int rows, int cols);
    /// Converts Matrix reference to SimpleMatrix
    SimpleMatrix(const Matrix& copy);
    /// Converts Matrix pointer to SimpleMatrix
    SimpleMatrix(const Matrix* copy);

    /// Destructor, frees memory
    ~SimpleMatrix();

    /// Creates an exact copy of the matrix and returns it.
    SimpleMatrix* clone() const;
    /// Copies cp's data onto this
    void copy(SimpleMatrix* cp);

    /// Set every element of this to val
    void set(double val);
    /// Copies lower triangle tri to matrix_
    void set(const double *tri);
    /// Set a single element of matrix_
    void set(int m, int n, double val) { matrix_[m][n] = val; }
    void set(SimpleVector *vec);
    void set(double **mat);

    /// Sets the diagonal of matrix_ to vec
    double get(int m, int n) { return matrix_[m][n]; }
    /// Returns matrix_
    double **to_block_matrix() const;
    /// Sets the name of the matrix
    void set_name(std::string name) { name_ = name; }

    /// Prints the matrix with print_mat
    void print(FILE *out = outfile);
    /// Print the matrix with corresponding eigenvalues below each column
    void eivprint(SimpleVector *values, FILE *out = outfile);
    /// The number of rows
    int rows() const { return rows_; }
    /// The number of columns
    int cols() const { return cols_; }
    /// Set matrix to identity
    void set_to_identity();
    /// Zero out the matrix
    void zero();
    /// Zero out the diagonal
    void zero_diagonal();

    /// Returns the trace of this
    double trace() const;
    /// Create a new SimpleMatrix which is the transpose of this
    SimpleMatrix *transpose();

    /// Add a matrix to this
    void add(const SimpleMatrix*);
    /// Subtracts a matrix from this
    void subtract(const SimpleMatrix*);
    /// Multiples the two arguments and adds their result to this
    void accumulate_product(const SimpleMatrix*, const SimpleMatrix*);
    /// Scales this matrix
    void scale(double);
    /// Returns the sum of the squares of this
    double sum_of_squares();
    /// Add val to an element of this
    void add(int m, int n, double val) { matrix_[m][n] += val; }
    /// Scale row m by a
    void scale_row(int m, double a);
    /// Scale column n by a
    void scale_column(int n, double a);
    /// Transform a by transformer save result to this
    void transform(SimpleMatrix* a, SimpleMatrix* transformer);
    /// Transform this by transformer
    void transform(SimpleMatrix* transformer);
    /// Back transform a by transformer save result to this
    void back_transform(SimpleMatrix* a, SimpleMatrix* transformer);
    /// Back transform this by transformer
    void back_transform(SimpleMatrix* transformer);

    /// Return the vector dot product of rhs by this
    double vector_dot(SimpleMatrix* rhs);

    /// General matrix multiply, saves result to this
    void gemm(bool transa, bool transb, double alpha, const SimpleMatrix* a, const SimpleMatrix* b, double beta);
    /// Diagonalize this, eigvector and eigvalues must be created by caller.
    void diagonalize(SimpleMatrix* eigvectors, SimpleVector* eigvalues);

    /// Saves the block matrix to PSIO object with fileno and with the toc position of the name of the matrix
    void save(Ref<psi::PSIO>& psio, unsigned int fileno);
    /// Saves the matrix in ASCII format to filename
    void save(const char *filename, bool append=true, bool saveLowerTriangle = true);
    /// Saves the matrix in ASCII format to filename
    void save(std::string filename, bool append=true, bool saveLowerTriangle = true);
};

//! Simple matrix reference wrapped class.
class RefSimpleMatrix : public Ref< SimpleMatrix > {
public:
    RefSimpleMatrix();
    RefSimpleMatrix(SimpleMatrix *o);
    RefSimpleMatrix(const RefSimpleMatrix& o);

    SimpleMatrix* clone() const;
    void copy(const RefSimpleMatrix& cp);

    void set(double val);
    void set(const double *tri);
    void set(int m, int n, double val);
    void set(const RefSimpleVector& vec);
    void set(double **mat);
    double get(int m, int n);
    double **to_block_matrix();

    void set_name(std::string name);

    void print(FILE *out = outfile);
    void eivprint(RefSimpleVector& values, FILE *out=outfile);
    int rows() const;
    int cols() const;
    void set_to_identity();
    void zero();
    void zero_diagonal();

    // Math routines
    double trace();
    RefSimpleMatrix transpose();

    void add(const RefSimpleMatrix&);
    void subtract(const RefSimpleMatrix&);
    void scale(double);
    double sum_of_squares();
    void add(int m, int n, double val);
    void scale_row(int m, double a);
    void scale_column(int n, double a);
    void transform(RefSimpleMatrix& a, RefSimpleMatrix& transformer);
    void transform(RefSimpleMatrix& transformer);
    void back_transform(RefSimpleMatrix& a, RefSimpleMatrix& transformer);
    void back_transform(RefSimpleMatrix& transformer);

    double vector_dot(RefSimpleMatrix& rhs);

    void gemm(bool transa, bool transb, double alpha, const RefSimpleMatrix& a, const RefSimpleMatrix& b, double beta);
    void diagonalize(RefSimpleMatrix& eigvectors, RefSimpleVector& eigvalues);

    void save(Ref<psi::PSIO>& psio, unsigned int fileno);
    void save(const char *filename, bool append=true, bool saveLowerTriangle = true);
    void save(std::string filename, bool append=true, bool saveLowerTriangle = true);

    // Make this refer to m
    RefSimpleMatrix& operator=(SimpleMatrix* m);
    // Make this and m refer to the same matrix
    RefSimpleMatrix& operator=(const RefSimpleMatrix& m);

    // Multiply this by a matrix and return a matrix
    RefSimpleMatrix operator*(const RefSimpleMatrix&) const;
    // Multiply this by a scalar and return the result
    RefSimpleMatrix operator*(double) const;
    // Matrix addition
    RefSimpleMatrix operator+(const RefSimpleMatrix&) const;
    // Matrix subtraction
    RefSimpleMatrix operator-(const RefSimpleMatrix&) const;
};

typedef Ref<RefMatrix, SimpleReferenceCount, StandardArrayPolicy> RefMatrixArray;
typedef Ref<RefVector, SimpleReferenceCount, StandardArrayPolicy> RefVectorArray;
typedef Ref<RefSimpleMatrix, SimpleReferenceCount, StandardArrayPolicy> RefSimpleMatrixArray;

#include "matrix_i.cc"

}

#endif // MATRIX_H
