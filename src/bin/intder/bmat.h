#ifndef __bmat_h
#define __bmat_h

class BMat
{
 
public:
  double **BMatrix;
  double **BMatSave;
  double *SVectArray;
  double **AMatrix;
  int disp;

  BMat();
  ~BMat();

  void init();

  void make_BMat();
  void invert_BMat();
  void BRow(double*, double, int, double*);
  void StoreElement(double*, int, int, double*, double*);
  void StoreElement(double*, int, int, int, double*, double*, double*);
  void StoreElement(double*, int, int, int, int, double*, double*, double*);
  void StoreElement(double*, int, int, int, int, double*, double*, double*, double*);
  void print_BMat();
};

#endif
