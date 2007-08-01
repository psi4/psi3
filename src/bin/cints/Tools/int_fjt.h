/*! \file int_fjt.h
    \ingroup (CINTS)
    \brief Enter brief description of file here 
*/

void init_fjt(int max);
void free_fjt();
void init_fjt_table(double_array_t *table);
void free_fjt_table(double_array_t *table);
void int_fjt(double_array_t *table, int J, double wval);

