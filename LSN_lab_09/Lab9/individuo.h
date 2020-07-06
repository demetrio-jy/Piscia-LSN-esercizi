#ifndef _individual_
#define _individual_

class Individual {

protected:
	int Nelem = 32;
	int* order = new int[Nelem];
	double L;

public:
	//constructor
	Individual();
	//destructor
	~Individual();
	//metodhs
	int GetCity(int);
	void SetCity(int, int);
	int* GetCities(void);
	void SetCities(int*);
	double GetL(void);
	void SetL(double);
};

#endif

