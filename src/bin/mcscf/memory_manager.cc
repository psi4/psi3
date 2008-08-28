#include <map>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <liboptions/liboptions.h>
#include "memory_manager.h"

extern FILE *infile, *outfile;

using namespace std;

MemoryManager::MemoryManager(){
  CurrentAllocated = 0;
  MaximumAllocated = 0;
  MaximumAllowed   = size_t(1024)*size_t(1024)*size_t(options_get_int("MEMORY"));
  allocated_memory    = 0.0;
  total_memory        = double(options_get_int("MEMORY"));
  integral_strip_size = total_memory * 0.05;
}

MemoryManager::~MemoryManager()
{
}


void MemoryManager::RegisterMemory(void *mem, AllocationEntry& entry, size_t size)
{
  AllocationTable[mem] = entry;
  CurrentAllocated += size;
  if (CurrentAllocated > MaximumAllocated)
    MaximumAllocated = CurrentAllocated;
  if(options_get_int("DEBUG") > 9){
    fprintf(outfile, "\n  ==============================================================================");
    fprintf(outfile, "\n  MemoryManager Allocated   %12ld bytes (%8.1f Mb)",size,double(size)/1048576.0);
    fprintf(outfile, "\n  %-15s allocated   at %s:%d", entry.variableName.c_str(), entry.fileName.c_str(), entry.lineNumber);
    fprintf(outfile, "\n  Currently used            %12ld bytes (%8.1f Mb)",CurrentAllocated,
                 double(CurrentAllocated)/1048576.0);
    fprintf(outfile, "\n  ==============================================================================");
    fflush(outfile);
  }
}

void MemoryManager::UnregisterMemory(void *mem, size_t size, const char *fileName, size_t lineNumber)
{
  CurrentAllocated -= size;
  AllocationEntry& entry = AllocationTable[mem];
  if(options_get_int("DEBUG") > 9){
    fprintf(outfile, "\n  ==============================================================================");
    fprintf(outfile, "\n  MemoryManager Deallocated %12ld bytes (%8.1f Mb)",size,double(size)/1048576.0);
    fprintf(outfile, "\n  %-15s allocated   at %s:%d", entry.variableName.c_str(), entry.fileName.c_str(), entry.lineNumber);
    fprintf(outfile, "\n  %-15s deallocated at %s:%d", entry.variableName.c_str(), fileName, lineNumber);
    fprintf(outfile, "\n  Currently used            %12ld bytes (%8.1f Mb)",CurrentAllocated,
                 double(CurrentAllocated)/1048576.0);
    fprintf(outfile, "\n  ==============================================================================");
    fflush(outfile);
  }
  AllocationTable.erase(mem);

}

void MemoryManager::MemCheck(FILE *output)
{
  static bool alreadyChecked = false;
  
  fprintf(output, "\n\n");
  fprintf(output, "  ==============================================================================\n");
  fprintf(output, "  Memory Usage Report\n\n");
  fprintf(output, "  Maximum memory used: %8.1f Mb \n",double(MaximumAllocated)/1048576.0);
  fprintf(output, "  Number of objects still in memory: %-6d  Current bytes used: %-12lu",CurrentAllocated,AllocationTable.size());

  if (AllocationTable.size() > 0) {
    if (alreadyChecked == false)
      fprintf(output, "\n\n  Attempting to free the following objects:\n");
    else
      fprintf(output, "\n\n  Unable to delete the following objects:\n");
                  
    std::map<void*, AllocationEntry>::iterator it;

    for (it=AllocationTable.begin(); it != AllocationTable.end(); it++){
      fprintf(output, "  %15s allocated at %s:%d\n", (*it).second.variableName.c_str(), (*it).second.fileName.c_str(), (*it).second.lineNumber);
      fflush(outfile);
    }
  }
  fprintf(output, "\n  ==============================================================================\n");
}

