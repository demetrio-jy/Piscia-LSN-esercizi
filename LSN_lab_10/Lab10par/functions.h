Random rnd;
int* seed = new int [4];
int p1;
int p2;

const int Ncities = 32;
const int Npopulation = 100;
double positions[Ncities][2];
double x[Ncities];
double y[Ncities];
int iterations = 10000;

int Nmigr = 1000;
int send[Ncities]{ 0 };
int receive[Ncities]{ 0 };

int method = 1; //0 = sulla circonferenza, 1 = dentro il quadrato 

Individual * population = new Individual[Npopulation];
Individual * nextGen = new Individual[Npopulation];
double* Lbest = new double[iterations + 1];
double* Lmean = new double[iterations + 1];
double* Lerr = new double[iterations + 1];
Individual best;
double Lmin = 100;
double prob;


void RandomInitializer(int);
void PositionGenerator(void);
void Coordinates(void);
void CoordInv(void);
int CheckFunction(Individual);
void PopulationGenerator(void);
double L(Individual);
void FitnessCalculator(void);
void Quicksort(Individual*, int, int);
Individual Selector(void);
void PairMutation(Individual);
void ShiftMutation(Individual);
void InversionMutation(Individual);
void PermutationMutation(Individual);
void CrossingOver(Individual, Individual);
int PBC(int);
double Mean(void);
double StdDev(void);
void CopySend(void);
void CopyReceive(void);
void Superfunction(int);