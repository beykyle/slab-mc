
void Input::setHist(int nHist) {
  nHistories = nHist;
}

void Input::setSlab(Slab slab_input) {
  slab = slab_input;
}

void Input::Input(void) {
  //create dictionary of inputs

  ifstream inpFile;
  inpFile.open("inp.i");
  if (inpFile.is_open())  {
    while ( getline (myfile,line) )  {
    	//split line at ":", strip whitespace on either side
        //add second as value of dictionary, with first element as key
    }
    myfile.close();
  }
