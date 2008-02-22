#ifndef _psi_src_bin_psimrcc_utilities_h
#define _psi_src_bin_psimrcc_utilities_h
/***************************************************************************
 *  PSIMRCC : Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/

#include <vector>
#include <iostream>
#include <sstream>
#include <sys/<ctime>>

namespace psi{ namespace psimrcc{

using namespace std;

typedef std::vector<std::string>            strvec;

bool space(char c);
bool not_space(char c);
std::vector<std::string> split(const std::string& str);
std::vector<std::string> split_indices(const std::string& str);
void to_lower(std::string& sBuffer);
std::string to_string(const int val);
std::string to_string(const double val);

//  cout << "\n\t" << k << "\t" << setiosflags(ios::fixed | ios::showpoint);
//                         cout << 
double StringTo(const std::string& s);
double ToDouble(const std::string inString);
int string_to_integer(const std::string inString);
void append_reference(std::string& str, int reference);
std::string add_reference(std::string& str, int reference);
void append_reference(std::string& str, int reference);

std::string find_and_replace(std::string & source, const std::string & target, const std::string & replace);
void trim_spaces(std::string& str);

/**
  @author Francesco Evangelista <frank@ccc.uga.edu>
*/
class Timer
{
public:
    Timer() {gettimeofday(&___start,&___dummy);}
    double get() {gettimeofday(&___end,&___dummy);
      delta_time_seconds=(___end.tv_sec - ___start.tv_sec) + (___end.tv_usec - ___start.tv_usec)/1000000.0;
      delta_time_hours=delta_time_seconds/3600.0;
      delta_time_days=delta_time_hours/24.0;
      return(delta_time_seconds);}
private:
  struct timeval ___start,___end;
  struct timezone ___dummy;
  double delta_time_seconds;
  double delta_time_hours;
  double delta_time_days;
};



//mise en place des variables, en début de fonction
#define SETUP_TIMING   struct timeval ___start,___end;struct timezone ___dummy;double delta_time_seconds;double delta_time_hours;double delta_time_days

//lancement du chrono
#define START_TIMING gettimeofday(&___start,&___dummy)

//arrêt du chrono et affichage
#define END_TIMING   gettimeofday(&___end,&___dummy); delta_time_seconds=(___end.tv_sec - ___start.tv_sec) + (___end.tv_usec - ___start.tv_usec)/1000000.0;


double to_MB(size_t n);

void print_dmatrix(double** matrix, int size1, int size2, FILE* outfile);
void print_dmatrix(double** matrix, int size1, int size2, FILE* outfile, char* text);
void print_dmatrix159(double** matrix, int size1, int size2, FILE* outfile, char* text);
void print_dmatrix106(double** matrix, int size1, int size2, FILE* outfile, char* text);

void scalar_multiply_dmatrix(double** A, double** B,double** C, int size1, int size2);
void scalar_divide_dmatrix(double** A, double** B,double** C, int size1, int size2);
void copy_dmatrix(double** A, double** B,int size1, int size2);
void zero_diagonal_dmatrix(double** A, int size1);
double dot_dmatrix(double** A, double** B, int size1, int size2);
int** init_2imatrix(int size1, int size2);
void free_2imatrix(int** matrix, int size1, int size2);

int*** init_3imatrix(int size1, int size2, int size3);
void free_3imatrix(int*** matrix, int size1, int size2, int size3);
double*** init_3dmatrix(int size1, int size2, int size3);
void free_3dmatrix(double*** matrix, int size1, int size2, int size3);
double ****init_4dmatrix(int size1, int size2, int size3, int size4);
void free_4dmatrix(double**** matrix, int size1, int size2, int size3, int size4);
void init_4dmatrix(double****& matrix,int size1, int size2, int size3, int size4);
void print_4dmatrix(double**** matrix, int size1, int size2, int size3, int size4, FILE* outfile, char* text);

unsigned long int init_smatrix(short**& matrix,int size1, int size2);
unsigned long int free_smatrix(short**& matrix, int size1, int size2);
unsigned long int init_smatrix(short***& matrix,int size1, int size2, int size3);
unsigned long int free_smatrix(short*** matrix, int size1, int size2, int size3);

unsigned long int init_imatrix(int**& matrix,int size1, int size2);
unsigned long int free_imatrix(int** matrix, int size1, int size2);
unsigned long int init_imatrix(int***& matrix,int size1, int size2, int size3);
unsigned long int free_imatrix(int*** matrix, int size1, int size2, int size3);

// template <typename T> size_t init_matrix(T*& matrix,size_t size);
// template <typename T> size_t free_matrix(T*& matrix,size_t size);
// template <typename T> size_t init_matrix(T**& matrix,size_t size1,size_t size2);
// template <typename T> size_t free_matrix(T**& matrix,size_t size1,size_t size2);
// template <typename T> size_t init_matrix(T***& matrix,size_t size1,size_t size2,size_t size3);
// template <typename T> size_t free_matrix(T***& matrix,size_t size1,size_t size2,size_t size3);
// 
// template <typename T>
// size_t init_matrix(T*& matrix,size_t size)
// {
//   if(size<=0){
//     matrix = NULL;
//   }else{
//     matrix    = new T[size];
//     for(size_t i=0;i<size;i++)
//       matrix[i]=static_cast<T>(0);   // Zero all the elements
//   }
//   return(size*sizeof(T));
// }
// 
// template <typename T>
// size_t free_matrix(T*& matrix,size_t size)
// {
//   if(matrix == NULL)
//     return(0);
//   delete[] matrix;
//   matrix = NULL;
//   return(size*sizeof(T));
// }
// 
// template <typename T>
// size_t init_matrix(T**& matrix,size_t size1,size_t size2)
// {
//   size_t size = size1*size2;
//   if(size<=0){
//     matrix = NULL;
//   }else{
//     matrix    = new T*[size1];
//     T* vector = new T[size];
//     for(size_t i=0;i<size;i++)
//       vector[i]=static_cast<T>(0);   // Zero all the elements
//     for(size_t i=0;i<size1;i++)
//       matrix[i]=&(vector[i*size2]);  // Assign the rows pointers
//   }
//   return(size*sizeof(T));
// }
// 
// template <typename T>
// size_t free_matrix(T**& matrix,size_t size1,size_t size2)
// {
//   size_t size = size1*size2;
//   if(matrix == NULL)
//     return(0);
//   delete[] matrix[0];
//   delete[] matrix;
//   matrix = NULL;
//   return(size*sizeof(T));
// }
// 
// template <typename T>
// size_t init_matrix(T***& matrix,size_t size1,size_t size2,size_t size3)
// {
//   size_t size = size1*size2*size3;
//   if(size<=0){
//     matrix = NULL;
//   }else{
//     matrix    = new T**[size1];
//     for(size_t i=0;i<size1;i++)
//       matrix[i]= new T*[size2];
//     T* vector = new T[size];
//     for(size_t i=0;i<size;i++)
//       vector[i]=static_cast<T>(0);   // Zero all the elements
//     for(size_t i=0;i<size1;i++)
//       for(size_t j=0;j<size2;j++)
//         matrix[i][j]=&(vector[i*size2*size3+j*size3]);  // Assign the rows pointers
//   }
//   return(size*sizeof(T));
// }
// 
// template <typename T>
// size_t free_matrix(T***& matrix,size_t size1,size_t size2,size_t size3)
// {
//   size_t size = size1*size2*size3;
//   if(matrix == NULL)
//     return(0);
//   delete[] matrix[0][0];
//   for(size_t i=0;i<size1;i++)
//     delete[] matrix[i];
//   delete[] matrix;
//   matrix = NULL;
//   return(size*sizeof(T));
// }

void generate_combinations(int n, int k, std::vector<std::vector<int> >& combinations);

}} /* End Namespaces */

#endif // _psi_src_bin_psimrcc_utilities_h