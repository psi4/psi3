//
// Reference Matrix class
//
inline RefMatrix::RefMatrix() {
    
}

inline RefMatrix::RefMatrix(Matrix *o) : Ref<Matrix>(o) {
    
}

inline RefMatrix::RefMatrix(const RefMatrix& o) : Ref<Matrix>(o) {
    
}

inline Matrix* RefMatrix::clone() const
{
    return new Matrix(this->pointer());
}

inline void RefMatrix::copy(const RefMatrix& cp)
{
    pointer()->copy(cp.pointer());
}

inline void RefMatrix::set(double val)
{
    pointer()->set(val);
}

inline void RefMatrix::set(const double *tri)
{
    pointer()->set(tri);
}

inline void RefMatrix::set(const double **sq)
{
    pointer()->set(sq);
}

inline void RefMatrix::set(int h, int m, int n, double val)
{
    pointer()->set(h, m, n, val);
}

inline void RefMatrix::set(const RefVector& vec)
{
    pointer()->set(vec.pointer());
}

inline void RefMatrix::set(const RefSimpleMatrix& sq)
{
    pointer()->set(sq.pointer());
}

inline double RefMatrix::get(int h, int m, int n)
{
    return pointer()->get(h, m, n);
}

inline double** RefMatrix::to_block_matrix()
{
    return pointer()->to_block_matrix();
}

inline void RefMatrix::set_name(std::string name)
{
    pointer()->set_name(name);
}

inline void RefMatrix::print(FILE *out, char *extra)
{
    pointer()->print(out, extra);
}

inline void RefMatrix::eivprint(RefVector& values, FILE *out)
{
    pointer()->eivprint(values.pointer(), out);
}

inline int* RefMatrix::rowspi() const
{
    return pointer()->rowspi();
}

inline int* RefMatrix::colspi() const
{
    return pointer()->colspi();
}

inline int RefMatrix::nirreps() const
{
    return pointer()->nirreps();
}

inline void RefMatrix::set_to_identity()
{
    pointer()->set_to_identity();
}

inline void RefMatrix::zero()
{
    pointer()->zero();
}

inline void RefMatrix::zero_diagonal()
{
    pointer()->zero_diagonal();
}

inline double RefMatrix::trace()
{
    return pointer()->trace();
}

inline RefMatrix RefMatrix::transpose()
{
    RefMatrix trans(pointer()->transpose());
    return trans;
}

inline RefSimpleMatrix RefMatrix::to_simple_matrix()
{
    RefSimpleMatrix simple(pointer()->to_simple_matrix());
    return simple;
}

inline void RefMatrix::add(const RefMatrix& rhs)
{
    pointer()->add(rhs.pointer());
}

inline void RefMatrix::subtract(const RefMatrix& rhs)
{
    pointer()->subtract(rhs.pointer());
}

inline void RefMatrix::scale(double val)
{
    pointer()->scale(val);
}

inline double RefMatrix::sum_of_squares()
{
    return pointer()->sum_of_squares();
}

inline void RefMatrix::add(int h, int m, int n, double val)
{
    pointer()->add(h, m, n, val);
}

inline void RefMatrix::scale_row(int h, int m, double a)
{
    pointer()->scale_row(h, m, a);
}

inline void RefMatrix::scale_column(int h, int n, double a)
{
    pointer()->scale_column(h, n, a);
}

inline void RefMatrix::transform(RefMatrix& a, RefMatrix& transformer)
{
    pointer()->transform(a.pointer(), transformer.pointer());
}

inline void RefMatrix::transform(RefMatrix& transformer)
{
    pointer()->transform(transformer.pointer());
}

inline void RefMatrix::back_transform(RefMatrix& a, RefMatrix& transformer)
{
    pointer()->back_transform(a.pointer(), transformer.pointer());
}

inline void RefMatrix::back_transform(RefMatrix& transformer)
{
    pointer()->back_transform(transformer.pointer());
}

inline double RefMatrix::vector_dot(RefMatrix& rhs)
{
    return pointer()->vector_dot(rhs.pointer());
}

inline void RefMatrix::gemm(bool transa, bool transb, double alpha, const RefMatrix& a, const RefMatrix& b, double beta)
{
    pointer()->gemm(transa, transb, alpha, a.pointer(), b.pointer(), beta);
}

inline void RefMatrix::diagonalize(RefMatrix& eigvectors, RefVector& eigvalues)
{
    pointer()->diagonalize(eigvectors.pointer(), eigvalues.pointer());
}

// Make this refer to m
inline RefMatrix& RefMatrix::operator=(Matrix* m)
{
    Ref<Matrix>::operator=(m);
    return *this;
}

// Make this and m refer to the same matrix
inline RefMatrix& RefMatrix::operator=(const RefMatrix& m)
{
    Ref<Matrix>::operator=(m);
    return *this;
}

inline RefMatrix RefMatrix::operator*(const RefMatrix&a) const
{
    RefMatrix r = a.clone();
    r->set(0.0);
    r->accumulate_product(pointer(), a.pointer());
    return r;
}

inline RefMatrix RefMatrix::operator*(double val) const
{
    RefMatrix r(clone());
    r.scale(val);
    return r;
}

inline RefMatrix RefMatrix::operator+(const RefMatrix&a) const
{
    RefMatrix ret(clone());
    ret->add(a.pointer());
    return ret;
}

inline RefMatrix RefMatrix::operator-(const RefMatrix&a) const
{
    RefMatrix ret(clone());
    ret->scale(-1.0);
    ret->add(a.pointer());
    return ret;
}

inline RefMatrix operator*(const double lhs, const RefMatrix& rhs)
{
    RefMatrix temp(rhs.clone());
    temp->scale(lhs);
    return temp;
}

inline bool RefMatrix::load(Ref<psi::PSIO>& psio, unsigned int fileno, char *tocentry, int nso)
{
    return pointer()->load(psio, fileno, tocentry, nso);
}

inline void RefMatrix::save(const char *filename, bool append, bool saveLowerTriangle, bool saveSubBlocks)
{
    pointer()->save(filename, append, saveLowerTriangle, saveSubBlocks);
}

inline void RefMatrix::save(std::string filename, bool append, bool saveLowerTriangle, bool saveSubBlocks)
{
    pointer()->save(filename, append, saveLowerTriangle, saveSubBlocks);
}

inline void RefMatrix::save(Ref<psi::PSIO>& psio, unsigned int fileno, bool saveSubBlocks)
{
    pointer()->save(psio, fileno, saveSubBlocks);
}

//
// Reference SimpleMatrix class
//
inline RefSimpleMatrix::RefSimpleMatrix() {
    
}

inline RefSimpleMatrix::RefSimpleMatrix(SimpleMatrix *o) : Ref<SimpleMatrix>(o) {
    
}

inline RefSimpleMatrix::RefSimpleMatrix(const RefSimpleMatrix& o) : Ref<SimpleMatrix>(o) {
    
}

inline SimpleMatrix* RefSimpleMatrix::clone() const
{
    return new SimpleMatrix(this->pointer());
}

inline void RefSimpleMatrix::copy(const RefSimpleMatrix& cp)
{
    pointer()->copy(cp.pointer());
}

inline void RefSimpleMatrix::set(double val)
{
    pointer()->set(val);
}

inline void RefSimpleMatrix::set(const double *tri)
{
    pointer()->set(tri);
}

inline void RefSimpleMatrix::set(int m, int n, double val)
{
    pointer()->set(m, n, val);
}

inline void RefSimpleMatrix::set(const RefSimpleVector& vec)
{
    pointer()->set(vec.pointer());
}

inline void RefSimpleMatrix::set(double **mat)
{
    pointer()->set(mat);
}

inline double RefSimpleMatrix::get(int m, int n)
{
    return pointer()->get(m, n);
}

inline double** RefSimpleMatrix::to_block_matrix()
{
    return pointer()->to_block_matrix();
}

inline void RefSimpleMatrix::set_name(std::string name)
{
    pointer()->set_name(name);
}

inline void RefSimpleMatrix::print(FILE *out)
{
    pointer()->print(out);
}

inline void RefSimpleMatrix::eivprint(RefSimpleVector& values, FILE *out)
{
    pointer()->eivprint(values.pointer(), out);
}

inline int RefSimpleMatrix::rows() const
{
    return pointer()->rows();
}

inline int RefSimpleMatrix::cols() const
{
    return pointer()->cols();
}

inline void RefSimpleMatrix::set_to_identity()
{
    pointer()->set_to_identity();
}

inline void RefSimpleMatrix::zero()
{
    pointer()->zero();
}

inline void RefSimpleMatrix::zero_diagonal()
{
    pointer()->zero_diagonal();
}

inline double RefSimpleMatrix::trace()
{
    return pointer()->trace();
}

inline RefSimpleMatrix RefSimpleMatrix::transpose()
{
    RefSimpleMatrix trans(pointer()->transpose());
    return trans;
}

inline void RefSimpleMatrix::add(const RefSimpleMatrix& rhs)
{
    pointer()->add(rhs.pointer());
}

inline void RefSimpleMatrix::subtract(const RefSimpleMatrix& rhs)
{
    pointer()->subtract(rhs.pointer());
}

inline void RefSimpleMatrix::scale(double val)
{
    pointer()->scale(val);
}

inline double RefSimpleMatrix::sum_of_squares()
{
    return pointer()->sum_of_squares();
}

inline void RefSimpleMatrix::add(int m, int n, double val)
{
    pointer()->add(m, n, val);
}

inline void RefSimpleMatrix::scale_row(int m, double a)
{
    pointer()->scale_row(m, a);
}

inline void RefSimpleMatrix::scale_column(int n, double a)
{
    pointer()->scale_column(n, a);
}

inline void RefSimpleMatrix::transform(RefSimpleMatrix& a, RefSimpleMatrix& transformer)
{
    pointer()->transform(a.pointer(), transformer.pointer());
}

inline void RefSimpleMatrix::transform(RefSimpleMatrix& transformer)
{
    pointer()->transform(transformer.pointer());
}

inline void RefSimpleMatrix::transform(const RefSimpleMatrix& transformer)
{
    pointer()->transform(transformer.pointer());
}

inline void RefSimpleMatrix::back_transform(RefSimpleMatrix& a, RefSimpleMatrix& transformer)
{
    pointer()->back_transform(a.pointer(), transformer.pointer());
}

inline void RefSimpleMatrix::back_transform(RefSimpleMatrix& transformer)
{
    pointer()->back_transform(transformer.pointer());
}

inline void RefSimpleMatrix::back_transform(const RefSimpleMatrix& transformer)
{
    pointer()->back_transform(transformer.pointer());
}

inline double RefSimpleMatrix::vector_dot(RefSimpleMatrix& rhs)
{
    return pointer()->vector_dot(rhs.pointer());
}

inline void RefSimpleMatrix::gemm(bool transa, bool transb, double alpha, const RefSimpleMatrix& a, const RefSimpleMatrix& b, double beta)
{
    pointer()->gemm(transa, transb, alpha, a.pointer(), b.pointer(), beta);
}

inline void RefSimpleMatrix::diagonalize(RefSimpleMatrix& eigvectors, RefSimpleVector& eigvalues)
{
    pointer()->diagonalize(eigvectors.pointer(), eigvalues.pointer());
}

inline void RefSimpleMatrix::save(Ref<psi::PSIO>& psio, unsigned int fileno)
{
    pointer()->save(psio, fileno);
}

inline void RefSimpleMatrix::save(const char *filename, bool append, bool saveLowerTriangle)
{
    pointer()->save(filename, append, saveLowerTriangle);
}

// Make this refer to m
inline RefSimpleMatrix& RefSimpleMatrix::operator=(SimpleMatrix* m)
{
    Ref<SimpleMatrix>::operator=(m);
    return *this;
}

// Make this and m refer to the same matrix
inline RefSimpleMatrix& RefSimpleMatrix::operator=(const RefSimpleMatrix& m)
{
    Ref<SimpleMatrix>::operator=(m);
    return *this;
}

inline RefSimpleMatrix RefSimpleMatrix::operator*(const RefSimpleMatrix&a) const
{
    RefSimpleMatrix r = a.clone();
    r->set(0.0);
    r->accumulate_product(pointer(), a.pointer());
    return r;
}

inline RefSimpleMatrix RefSimpleMatrix::operator*(double val) const
{
    RefSimpleMatrix r(clone());
    r.scale(val);
    return r;
}

inline RefSimpleMatrix RefSimpleMatrix::operator+(const RefSimpleMatrix&a) const
{
    RefSimpleMatrix ret(clone());
    ret->add(a.pointer());
    return ret;
}

inline RefSimpleMatrix RefSimpleMatrix::operator-(const RefSimpleMatrix&a) const
{
    RefSimpleMatrix ret(clone());
    ret->scale(-1.0);
    ret->add(a.pointer());
    return ret;
}

inline RefSimpleMatrix operator*(const double lhs, const RefSimpleMatrix& rhs)
{
    RefSimpleMatrix temp(rhs.clone());
    temp->scale(lhs);
    return temp;
}
