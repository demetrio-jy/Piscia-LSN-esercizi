Random rnd;
int* seed = new int [4];
int p1;
int p2;

const int Ncities = 32;
const int Npopulation = 200;
double positions[Ncities][2];
int iterations = 20000;

int method = 1; //0 = sulla circonferenza, 1 = dentro il quadrato 

Individual * population = new Individual[Npopulation];
double* Lbest = new double[iterations + 1];
Individual best;
double Lmin = 100;
double prob;

double beta; //temperature for simulated annealing
double Tmax = 50; //Tmin = 1
int index;


void RandomInitializer(void);
void PositionGenerator(void);
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
void Superfunction(int);