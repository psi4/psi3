#ifndef _psi_src_lib_libutil_libutil_h_
#define _psi_src_lib_libutil_libutil_h_

#include <vector>
#include <iostream>
#include <sstream>
#include <sys/time.h>

extern FILE *outfile;

using namespace std;

typedef std::vector<std::string>            strvec;

bool space(char c);
bool not_space(char c);
std::vector<std::string> split(const std::string& str);
std::vector<std::string> split_indices(const std::string& str);
void to_lower(std::string& sBuffer);
std::string to_string(const int val);
std::string to_string(const double val);

double StringTo(const std::string& s);
double ToDouble(const std::string inString);
int string_to_integer(const std::string inString);
void append_reference(std::string& str, int reference);
std::string add_reference(std::string& str, int reference);
void append_reference(std::string& str, int reference);

std::string find_and_replace(std::string & source, const std::string & target, const std::string & replace);
void trim_spaces(std::string& str);

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

double to_MB(size_t n);

void print_dmatrix(double** matrix, int size1, int size2, FILE* outfile);
void print_dmatrix(double** matrix, int size1, int size2, FILE* outfile, const char* text);
void print_dmatrix159(double** matrix, int size1, int size2, FILE* outfile, const char* text);
void print_dmatrix106(double** matrix, int size1, int size2, FILE* outfile, const char* text);

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

void print_error(std::string message, const char* file, int line);
void print_error(const char* message, const char* file, int line);
void print_error(const char* message, const char* file, int line,int error);
void print_developing(const char* message, const char* file, int line);
void print_developing(const char* message, const char* file, int line,int error);

void generate_combinations(int n, int k, std::vector<std::vector<int> >& combinations);

#endif // _psi_src_lib_libutil_libutil_h_

