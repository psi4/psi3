#ifndef _psi_src_bin_psimrcc_calculation_options_h
#define _psi_src_bin_psimrcc_calculation_options_h
/***************************************************************************
 *  PSIMRCC
 *  Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/

#include <map>
#include <string>

namespace psi{ namespace psimrcc{

typedef struct 
{
  bool option;
} CalculationBoolOption;

typedef struct 
{
  int option;
} CalculationIntOption;

typedef struct 
{
  double option;
} CalculationDoubleOption;

typedef struct 
{
  std::string option;
  std::string choices;
} CalculationStringOption;

typedef std::map<std::string,CalculationBoolOption>       BoolOptionsMap;
typedef std::map<std::string,CalculationIntOption>         IntOptionsMap;
typedef std::map<std::string,CalculationDoubleOption>   DoubleOptionsMap;
typedef std::map<std::string,CalculationStringOption>   StringOptionsMap;

class CalculationOptions
{
public:
  CalculationOptions();
  ~CalculationOptions();

  void        print();
  void        read_options();

  void        add_bool_option(char* cstr_option,bool bool_default);
  void        add_int_option(char* cstr_option,int int_default);
  void        add_double_option(char* cstr_option,double double_default);
  void        add_str_option(char* cstr_option,char* cstr_default);
  void        add_str_option_with_choices(char* cstr_option,char* cstr_default,char* cstr_choices);

  bool        get_bool_option(char* cstr_option);
  int         get_int_option(char* cstr_option);
  double      get_double_option(char* cstr_option);
  std::string get_str_option(char* cstr_option);

  void        set_bool_option(char* cstr_option, bool bool_value);
  void        set_int_option(char* cstr_option,int int_value);
  void        set_double_option(char* cstr_option,double double_value);
  void        set_str_option(char* cstr_option,char* cstr_value);
private:
  void read_bool(BoolOptionsMap::iterator& it);
  void read_int(IntOptionsMap::iterator& it);
  void read_double(DoubleOptionsMap::iterator& it);
  void read_string(StringOptionsMap::iterator& it);


  BoolOptionsMap      bool_options;
  IntOptionsMap       int_options;
  DoubleOptionsMap double_options;
  StringOptionsMap string_options;
};

extern CalculationOptions *options;

#endif

}} /* End Namespaces */