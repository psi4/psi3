/*#################################################################################
  coord_base.h

  coordinate base class declaration and definition
  ###############################################################################*/

template<class coord_type>class coord_base {
  public:
    coord_type* coord_arr;
    coord_base( int num_coord ) {
      coord_arr = new coord_type[num_coord];
    }
    ~coord_base() { delete [] coord_arr; }
};

