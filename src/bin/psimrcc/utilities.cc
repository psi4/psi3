/***************************************************************************
 *  PSIMRCC : Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/

#include "utilities.h"
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <algorithm>

namespace psi{ namespace psimrcc{



/***************************************************************************
 String manipulation
 ***************************************************************************/

std::vector<std::string> split(const std::string& str){
  // Split a string
  typedef string::const_iterator iter;
  strvec splitted_string;
  iter i = str.begin();
  while(i != str.end()){
    // Ignore leading blanks
    i = find_if(i,str.end(), not_space);
    // Find the end of next word
    iter j = find_if(i,str.end(),space);
    // Copy the characters in [i,j)
    if(i!=str.end())
      splitted_string.push_back(string(i,j));
    i = j;
  }
  return(splitted_string);
}


bool opening_square_bracket(char c);
bool closing_square_bracket(char c);

std::vector<std::string> split_indices(const std::string& str){
  // Split a string
  typedef string::const_iterator iter;
  strvec splitted_string;
  iter i = str.begin();
  while(i != str.end()){
    // Ignore leading blanks
    i = find_if(i,str.end(), opening_square_bracket);
    // Find the end of next word
    iter j = find_if(i,str.end(),closing_square_bracket);
    // Copy the characters in [i,j]
    if(i!=str.end())
      splitted_string.push_back(string(i,j+1));
    i = j;
  }
  return(splitted_string);
}

bool opening_square_bracket(char c)
{
  return (c=='[');
}

bool closing_square_bracket(char c)
{
  return (c==']');
}

bool space(char c)
{
  return isspace(c);
}

bool not_space(char c)
{
  return !isspace(c);
}

void to_lower(std::string& sBuffer)
{
  std::transform( sBuffer.begin(), sBuffer.end(), sBuffer.begin(),::tolower);
}

double ToDouble(const std::string inString)
{
  double d = 0;
  char* end;
  d = std::strtod( inString.c_str(), &end ); 
  return d;
}

int string_to_integer(const std::string inString)
{
  int i = 0;
  char* end;
  i = static_cast<int>(std::strtod( inString.c_str(), &end )); 
  return i;
}

std::string add_reference(std::string& str, int reference)
{
  return(str + "{" + to_string(reference) + "}");
}

void append_reference(std::string& str, int reference)
{
  str += "{" + to_string(reference) + "}";
}

std::string to_string(const int val)
{
    std::stringstream strm;
    strm <<  val;
    return strm.str();
}

std::string to_string(const double val)
{
    std::stringstream strm;
    strm << setprecision(25) << setw(35)  << val;
    return strm.str();
}

double StringTo(const std::string& s)
{
    std::stringstream ss;
    ss << s;
    double x;
    ss >> x;
    return x;
}

std::string find_and_replace(std::string & source, const std::string & target, const std::string & replace )
{
  string str = source;
  string::size_type pos = 0;   // where we are now
  string::size_type found;     // where the found data is

  if (target.size () > 0)   // searching for nothing will cause a loop
    {
    while ((found = str.find (target, pos)) != string::npos)
      {
      str.replace (found, target.size (), replace);
      pos = found + replace.size ();
      }
    }
  return str;
}

void trim_spaces(std::string& str)
{
  // Trim Both leading and trailing spaces
  size_t startpos = str.find_first_not_of(" \t"); // Find the first character position after excluding leading blank spaces
  size_t endpos = str.find_last_not_of(" \t"); // Find the first character position from reverse af

  // if all spaces or empty return an empty string
  if(( string::npos == startpos ) || ( string::npos == endpos))
  {
    str = "";
  }
  else
    str = str.substr( startpos, endpos-startpos+1 );
}

/*********************************************************
  Combinatorial Functions
 *********************************************************/

void generate_combinations(int n, int k, std::vector<std::vector<int> >& combinations)
{
  if((n>0)&&(k>0)){
    int i,t,index;
    int *a,*comb,*list;
    std::vector<int> combination;

    list  = new int[n];
    a     = new int[k];
    comb  = new int[k];

    for(i=0;i<n;i++)
      list[i]=i;

    index=0;
    t=0;
    a[0]=0;
    loop3: comb[t]=list[a[t]];
    if(t==k-1)
      goto loop12;
    a[t+1]=a[t]+1;
    if (a[t+1]==n)
      goto loop9;
    t++;
    goto loop3;

    loop9: t--;
    if (t==-1)
      goto loop16;
    goto loop13;

    loop12: index++;

    // Store the combination
    combination.clear();
    for(i=0;i<k;i++)
      combination.push_back(comb[i]);
    combinations.push_back(combination);

    loop13: a[t]++;
    if (a[t]==n)
      goto loop9;
    goto loop3;
    loop16: delete[] a;
    delete[] comb;
    delete[] list;
  }
}


/*********************************************************
 Memory Allocation
 *********************************************************/

/**
 * Convert the size of a doubles array in Mb using the definition 1Mb = 1048576 bytes
 * @param n size of the array
 * @return 
 */
double to_MB(size_t n)
{
  return(double(n*sizeof(double))/1048576.0);
  // Using this definition 1 Mb has ca. 5% more than 1000000 bytes
}

void scalar_multiply_dmatrix(double** A, double** B,double** C, int size1, int size2)
{
  int i,j;
  double *a,*b,*c; 
  for(i=0;i<size1;i++){
    a = A[i];
    b = B[i];
    c = C[i];
    for(j=0;j<size2;j++)
      c[j]=a[j]*b[j];
  }
}


void scalar_divide_dmatrix(double** A, double** B,double** C, int size1, int size2)
{
  int i,j;
  double *a,*b,*c; 
  for(i=0;i<size1;i++){
    a = A[i];
    b = B[i];
    c = C[i];
    for(j=0;j<size2;j++)
      c[j]=a[j]/b[j];
  }
}


void copy_dmatrix(double** A, double** B,int size1, int size2)
{
  // Copy A in B
  int i,j;
  double *a,*b; 
  for(i=0;i<size1;i++){
    a = A[i];
    b = B[i];
    for(j=0;j<size2;j++)
      b[j]=a[j];
  }
}


void zero_diagonal_dmatrix(double** A, int size1)
{
  // Zero A[i][i]
  int i;
  for(i=0;i<size1;i++)
    A[i][i]=0.0;
}


double dot_dmatrix(double** A, double** B, int size1, int size2)
{
  int i,j;
  double *a,*b;
  double sum=0.0; 
  for(i=0;i<size1;i++){
    a = A[i];
    b = B[i];
    for(j=0;j<size2;j++)
      sum+=a[j]*b[j];
  }
  return(sum);
}

void print_dmatrix(double** matrix, int size1, int size2, FILE* outfile)
{
  for(int i=0;i<size1;i++){
    fprintf(outfile,"\n");
    for(int j=0;j<size2;j++)
        fprintf(outfile,"%20.15f ",matrix[i][j]);
  }
  fprintf(outfile,"\n");
}

void print_dmatrix(double** matrix, int size1, int size2, FILE* outfile, char* text)
{
    fprintf(outfile,"\n\n\t%s\n",text);
    for(int i=0;i<size1;i++){
      fprintf(outfile,"\n");
      for(int j=0;j<size2;j++)
         fprintf(outfile,"%20.15f ",matrix[i][j]);
    }
  fprintf(outfile,"\n");
}

void print_dmatrix159(double** matrix, int size1, int size2, FILE* outfile, char* text)
{
    fprintf(outfile,"\n\n\t%s\n",text);
    for(int i=0;i<size1;i++){
      fprintf(outfile,"\n\t");
      for(int j=0;j<size2;j++)
         fprintf(outfile,"%15.9f",matrix[i][j]);
    }
  fprintf(outfile,"\n");
}

void print_dmatrix106(double** matrix, int size1, int size2, FILE* outfile, char* text)
{
    fprintf(outfile,"\n\n\t%s\n",text);
    for(int i=0;i<size1;i++){
      fprintf(outfile,"\n\t");
      for(int j=0;j<size2;j++)
         fprintf(outfile,"%10.6f",matrix[i][j]);
    }
  fprintf(outfile,"\n");
}

int** init_2imatrix(int size1, int size2)
{
    int** matrix = new int*[size1];
    for(int i=0;i<size1;i++){
        matrix[i] = new int[size2];
        for(int j=0;j<size2;j++)
            matrix[i][j]=0;
    }
    return(matrix);
}

void free_2imatrix(int** matrix, int size1, int size2)
{
    for(int i=0;i<size1;i++){
        delete[] matrix[i];
    }
    delete[] matrix;
}

int*** init_3imatrix(int size1, int size2, int size3)
{
    int*** matrix = new int**[size1];
    for(int i=0;i<size1;i++){
        matrix[i] = new int*[size2];
    }
    for(int i=0;i<size1;i++){
        for(int j=0;j<size2;j++){
            matrix[i][j] = new int[size3];
        }
    }
    return(matrix);
}

void free_3imatrix(int*** matrix, int size1, int size2, int size3)
{
    for(int i=0;i<size1;i++){
        for(int j=0;j<size2;j++){
            delete[] matrix[i][j];
        }
    }
    for(int i=0;i<size1;i++){
        delete[] matrix[i];
    }
    delete[] matrix;
}

double*** init_3dmatrix(int size1, int size2, int size3)
{
    double*** matrix;
    if((size1 ==0)||(size2 ==0)||(size3 ==0)){
      matrix=NULL;
    }else{
      matrix = new double**[size1];
      for(int i=0;i<size1;i++){
          matrix[i] = new double*[size2];
      }
      for(int i=0;i<size1;i++){
          for(int j=0;j<size2;j++){
              matrix[i][j] = new double[size3];
              for(int k=0;k<size3;k++){
                  matrix[i][j][k]=0.0;
              }
          }
      }
    }
    return(matrix);
}

void free_3dmatrix(double*** matrix, int size1, int size2, int size3)
{
  if(matrix!=NULL){
    for(int i=0;i<size1;i++){
        for(int j=0;j<size2;j++){
            delete[] matrix[i][j];
        }
    }
    for(int i=0;i<size1;i++){
        delete[] matrix[i];
    }
    delete[] matrix;
  }
}

void init_4dmatrix(double****& matrix,int size1, int size2, int size3, int size4)
{
    if((size1 ==0)||(size2 ==0)||(size3 ==0)||(size4 ==0)){
      printf("\n\n\tNULL Matrix\n");
      matrix=NULL;
    }else{
    matrix = new double***[size1];
    for(int i=0;i<size1;i++){
        matrix[i] = new double**[size2];
    }
    for(int i=0;i<size1;i++){
        for(int j=0;j<size2;j++){
            matrix[i][j] = new double*[size3];
        }
    }
    for(int i=0;i<size1;i++){
        for(int j=0;j<size2;j++){
            for(int k=0;k<size3;k++){
                matrix[i][j][k] = new double[size4];
                for(int l=0;l<size4;l++){
                    matrix[i][j][k][l]=0.0;
                }
            }
        }
    }
  }
  //printf("\n\tAllocated %9d MB",size1*size2*size3*size4*sizeof(double)/1048576);
}

double ****init_4dmatrix(int size1, int size2, int size3, int size4)
{

    double**** matrix;
    if((size1 ==0)||(size2 ==0)||(size3 ==0)||(size4 ==0)){
      printf("\n\n\tNULL Matrix\n");
      matrix=NULL;
    }else{
    matrix = new double***[size1];
    for(int i=0;i<size1;i++){
        matrix[i] = new double**[size2];
    }
    for(int i=0;i<size1;i++){
        for(int j=0;j<size2;j++){
            matrix[i][j] = new double*[size3];
        }
    }
    for(int i=0;i<size1;i++){
        for(int j=0;j<size2;j++){
            for(int k=0;k<size3;k++){
                matrix[i][j][k] = new double[size4];
                for(int l=0;l<size4;l++){
                    matrix[i][j][k][l]=0.0;
                }
            }
        }
    }
  }
  // printf("\n\tAllocated %9d MB",size1*size2*size3*size4*sizeof(double)/1048576);
    return(matrix);
}

void free_4dmatrix(double**** matrix, int size1, int size2, int size3, int size4)
{
  if((size1 ==0)||(size2 ==0)||(size3 ==0)||(size4 ==0)){
      printf("\n\n\tNULL Matrix\n");
  }
  if(matrix!=NULL){
    for(int i=0;i<size1;i++){
        for(int j=0;j<size2;j++){
            for(int k=0;k<size3;k++){
                delete[] matrix[i][j][k];
            }
        }
    }
    for(int i=0;i<size1;i++){
        for(int j=0;j<size2;j++){
                delete[] matrix[i][j];
        }
    }
    for(int i=0;i<size1;i++){
        delete[] matrix[i];
    }
    delete[] matrix;
    }
}

void print_4dmatrix(double**** matrix, int size1, int size2, int size3, int size4, FILE* outfile, char* text)
{
    fprintf(outfile,"\n\n\t%s\n",text);
    for(int i=0;i<size1;i++)
      for(int j=0;j<size2;j++)
        for(int k=0;k<size3;k++)
          for(int l=0;l<size4;l++)
            if(fabs(matrix[i][j][k][l])>1.0e-8)
              fprintf(outfile,"\n\t[%3d,%3d,%3d,%3d] -> %20.15f ",i,j,k,l,matrix[i][j][k][l]);
}

unsigned long int init_smatrix(short**& matrix,int size1, int size2)
{
  unsigned long int size,uli_size1,uli_size2;
  uli_size1 = (unsigned long int)size1;
  uli_size2 = (unsigned long int)size2;
  size=uli_size1*uli_size2;
  if(!uli_size1 || !uli_size2){
    matrix=NULL;
  }else{
    matrix = new short*[uli_size1];
    short* vector = new short[size];
    for(unsigned long int i=0;i<size;i++) vector[i]=0;
    for(unsigned long int i=0;i<uli_size1;i++)
      matrix[i]=&(vector[i*uli_size2]);
  }
  return(size*sizeof(short));
}

unsigned long int free_smatrix(short**& matrix, int size1, int size2)
{
  unsigned long int size,uli_size1,uli_size2;
  uli_size1 = (unsigned long int)size1;
  uli_size2 = (unsigned long int)size2;
  size=uli_size1*uli_size2;
  if(matrix == NULL) return(0);
  delete[] matrix[0];
  delete[] matrix;
  return(size*sizeof(short));
}

unsigned long int init_smatrix(short***& matrix,int size1, int size2, int size3)
{
  unsigned long int size = (unsigned long int)size1*size2*size3;
  matrix = new short**[size1];
  for(int i=0;i<size1;i++){
      matrix[i] = new short*[size2];
  }
  for(int i=0;i<size1;i++){
      for(int j=0;j<size2;j++){
          matrix[i][j] = new short[size3];
      }
  }
  return(size*sizeof(short));
}

unsigned long int free_smatrix(short*** matrix, int size1, int size2, int size3)
{
  unsigned long int size = (unsigned long int)size1*size2*size3;
  for(int i=0;i<size1;i++){
      for(int j=0;j<size2;j++){
          delete[] matrix[i][j];
      }
  }
  for(int i=0;i<size1;i++){
      delete[] matrix[i];
  }
  delete[] matrix;
  return(size*sizeof(short));
}


unsigned long int init_imatrix(int**& matrix,int size1, int size2)
{
  unsigned long int size,uli_size1,uli_size2;
  uli_size1 = (unsigned long int)size1;
  uli_size2 = (unsigned long int)size2;
  size=uli_size1*uli_size2;
  if(!uli_size1 || !uli_size2){
    matrix=(int **) NULL;
  }else{
    matrix = new int*[uli_size1];
    int*  vector = new int[size];
    for(unsigned long int i=0;i<size;i++) vector[i]=0;
    for(unsigned long int i=0;i<uli_size1;i++)
      matrix[i]=&(vector[i*uli_size2]);
  }
  return(size*sizeof(int));
}

unsigned long int free_imatrix(int** matrix, int size1, int size2)
{
  unsigned long int size,uli_size1,uli_size2;
  uli_size1 = (unsigned long int)size1;
  uli_size2 = (unsigned long int)size2;
  size=uli_size1*uli_size2;
  if(matrix == NULL) return(0);
  delete[] matrix[0];
  delete[] matrix;
  return(size*sizeof(int));
}

unsigned long int init_imatrix(int***& matrix,int size1, int size2, int size3)
{
  unsigned long int size = (unsigned long int)size1*size2*size3;
  matrix = new int**[size1];
  for(int i=0;i<size1;i++){
      matrix[i] = new int*[size2];
  }
  for(int i=0;i<size1;i++){
      for(int j=0;j<size2;j++){
          matrix[i][j] = new int[size3];
      }
  }
  return(size*sizeof(int));
}

unsigned long int free_imatrix(int*** matrix, int size1, int size2, int size3)
{
  unsigned long int size = (unsigned long int)size1*size2*size3;
  for(int i=0;i<size1;i++){
      for(int j=0;j<size2;j++){
          delete[] matrix[i][j];
      }
  }
  for(int i=0;i<size1;i++){
      delete[] matrix[i];
  }
  delete[] matrix;
  return(size*sizeof(int));
}

}} /* End Namespaces */