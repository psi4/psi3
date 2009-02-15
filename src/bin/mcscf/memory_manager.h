#ifndef _psi_src_bin_psimrcc_memory_manager_h
#define _psi_src_bin_psimrcc_memory_manager_h

#include <map>
#include <vector>
#include <string>

namespace psi{ namespace mcscf{

typedef struct {
	void *variable;
	std::string type;
	std::string variableName;
	std::string fileName;
	size_t lineNumber;
	std::vector<size_t> argumentList;
} AllocationEntry;

class MemoryManager
{
public:
  MemoryManager();
  ~MemoryManager();

  void MemCheck(FILE *output);

  double get_total_memory() const {return(total_memory);}
  // Memory handling routines
  void        add_allocated_memory(double value)       {allocated_memory+=value;}
  double      get_free_memory()                  const {return(total_memory-allocated_memory);}
  double      get_allocated_memory()             const {return(allocated_memory);}
  double      get_integral_strip_size()          const {return(integral_strip_size);}
  size_t      get_FreeMemory()                   const {return(MaximumAllowed - CurrentAllocated);}

  template <typename T>
  void allocate(const char *type, T*& matrix, size_t size, const char *variableName, const char *fileName, size_t lineNumber);
  template <typename T>
  void release_one(T*& matrix, const char *fileName, size_t lineNumber);

  template <typename T>
  void allocate(const char *type, T**& matrix, size_t size1, size_t size2, const char *variableName, const char *fileName, size_t lineNumber);
  template <typename T>
  void release_two(T**& matrix, const char *fileName, size_t lineNumber);

  template <typename T>
  void allocate(const char *type, T***& matrix,size_t size1,size_t size2,size_t size3, const char *variableName, const char *fileName, size_t lineNumber);
  template <typename T>
  void release_three(T***& matrix, const char *fileName, size_t lineNumber);
private:
  void RegisterMemory(void *mem, AllocationEntry& entry, size_t size);
  void UnregisterMemory(void *mem, size_t size, const char *fileName, size_t lineNumber);

  size_t CurrentAllocated;
  size_t MaximumAllocated;
  size_t MaximumAllowed;
  std::map<void *, AllocationEntry> AllocationTable;

  double total_memory;
  double allocated_memory;
  double integral_strip_size;
};

extern MemoryManager *mem;

template <typename T>
void MemoryManager::allocate(const char *type, T*& matrix, size_t size, const char *variableName, const char *fileName, size_t lineNumber)
{
  AllocationEntry newEntry;
	
  if(size<=0){
    matrix = NULL;
  }else{
    matrix    = new T[size];
    for(size_t i=0;i<size;i++)
      matrix[i]=static_cast<T>(0);   // Zero all the elements

    newEntry.variable = matrix;
    newEntry.type = type;
    newEntry.variableName = variableName;
    newEntry.fileName = fileName;
    newEntry.lineNumber = lineNumber;
    newEntry.argumentList.push_back(size);
    RegisterMemory((void*)matrix, newEntry, size*sizeof(T));
  }
}

template <typename T>
void MemoryManager::release_one(T*& matrix, const char *fileName, size_t lineNumber)
{
  if(matrix == NULL)
    return;

  size_t size = AllocationTable[(void*)matrix].argumentList[0];
  
  UnregisterMemory((void*)matrix, size*sizeof(T),fileName,lineNumber);
	
  delete[] matrix;
  matrix = NULL;
}

template <typename T>
void MemoryManager::allocate(const char *type, T**& matrix, size_t size1, size_t size2, const char *variableName, const char *fileName, size_t lineNumber)
{
  AllocationEntry newEntry;
  size_t size = size1*size2;

  if(size<=0){
    matrix = NULL;
    return;
  }else{
    matrix    = new T*[size1];
    T* vector = new T[size];
    for(size_t i=0;i<size;i++)
      vector[i]=static_cast<T>(0);   // Zero all the elements
    for(size_t i=0;i<size1;i++)
      matrix[i]=&(vector[i*size2]);  // Assign the rows pointers

    newEntry.variable = matrix;
    newEntry.type = type;
    newEntry.variableName = variableName;
    newEntry.fileName = fileName;
    newEntry.lineNumber = lineNumber;
    newEntry.argumentList.push_back(size1);
    newEntry.argumentList.push_back(size2);
    RegisterMemory((void*)matrix, newEntry, size*sizeof(T));
  }
}

template <typename T>
void MemoryManager::release_two(T**& matrix, const char *fileName, size_t lineNumber)
{
  if(matrix == NULL)
    return;

  size_t size = AllocationTable[(void*)matrix].argumentList[0] * AllocationTable[(void*)matrix].argumentList[1];

  UnregisterMemory((void*)matrix, size*sizeof(T),fileName,lineNumber);

  delete[] matrix[0];
  delete[] matrix;
  matrix = NULL;
}

template <typename T>
void MemoryManager::allocate(const char *type, T***& matrix,size_t size1,size_t size2,size_t size3, const char *variableName, const char *fileName, size_t lineNumber)
{
  AllocationEntry newEntry;
  size_t size = size1*size2*size3;
  if(size<=0){
    matrix = NULL;
    return;
  }else{
    matrix    = new T**[size1];
    for(size_t i=0;i<size1;i++)
      matrix[i]= new T*[size2];
    T* vector = new T[size];
    for(size_t i=0;i<size;i++)
      vector[i]=static_cast<T>(0);   // Zero all the elements
    for(size_t i=0;i<size1;i++)
      for(size_t j=0;j<size2;j++)
        matrix[i][j]=&(vector[i*size2*size3+j*size3]);  // Assign the rows pointers

    newEntry.variable = matrix;
    newEntry.type = type;
    newEntry.variableName = variableName;
    newEntry.fileName = fileName;
    newEntry.lineNumber = lineNumber;
    newEntry.argumentList.push_back(size1);
    newEntry.argumentList.push_back(size2);
    newEntry.argumentList.push_back(size3);
    RegisterMemory((void*)matrix, newEntry, size*sizeof(T));
  }
}

template <typename T>
void MemoryManager::release_three(T***& matrix, const char *fileName, size_t lineNumber)
{
  if(matrix == NULL)
    return;

  size_t size1 = AllocationTable[(void*)matrix].argumentList[0];
  size_t size = AllocationTable[(void*)matrix].argumentList[0] * AllocationTable[(void*)matrix].argumentList[1]
              * AllocationTable[(void*)matrix].argumentList[2];

  UnregisterMemory((void*)matrix, size*sizeof(T),fileName,lineNumber);

  delete[] matrix[0][0];
  for(size_t i=0;i<size1;i++)
    delete[] matrix[i];
  delete[] matrix;
  matrix = NULL;
}

#define allocate1(type, variable, size) \
	mem->allocate(#type, variable, size, #variable, __FILE__, __LINE__);
#define release1(variable) \
	mem->release_one(variable, __FILE__, __LINE__);
	
#define allocate2(type, variable, size1, size2) \
	mem->allocate(#type, variable, size1, size2, #variable, __FILE__, __LINE__);
#define release2(variable) \
	mem->release_two(variable, __FILE__, __LINE__);

#define allocate3(type, variable, size1, size2, size3) \
	mem->allocate(#type, variable, size1, size2, size3, #variable, __FILE__, __LINE__);
#define release3(variable) \
	mem->release_three(variable, __FILE__, __LINE__);

}}

#endif // _psi_src_bin_psimrcc_memory_manager_h
