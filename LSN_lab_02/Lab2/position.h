#ifndef __Position__
#define __Position__

class Position {

private:
  double p_x, p_y, p_z;

protected:

public:
  // constructors
  Position();
  // destructor
  ~Position();
  // methods
  void SetX(double);
  void SetY(double);
  void SetZ(double);
  double GetX();
  double GetY();
  double GetZ(); 
  double GetR();
};

#endif // __Position__
