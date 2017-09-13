
class Input {
  private:  
    Slab    slab;
    double  nHistories;
  public:
    Input(std::string fname);
      //constructor that reads from input file
      //sets input parameters, creates slab object
      void setHist(int nHist);
      void Input(void);
      void setSlab(Slab slab);
};

