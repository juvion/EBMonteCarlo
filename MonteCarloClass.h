//Use the precompiler to make sure the class definition is not included more than once.
#if !defined(MONTECARLO_H)
#define MONTECARLO_H
//#define REPLICADEBUG
#include "../RNAstructureZero/RNA_class/RNA.h"
#include "../RNAstructureZero/src/random.h"
#include <vector>
#include <iomanip>

#ifdef REPLICADEBUG
#include <fstream>
#endif

//STAS
//Counts the number of iterations while the program is running
// #define ItCounter 
#undef ItCounter

//#define REPLICASTATS
#undef REPLICASTATS

//If "ListOfPairs" is defined a list of canonical pairs is generated. 
//#define ListOfPairs
//Equal to: Steps to form NonCanonical pairs are rejected and not counted. (NCan_RN)

//If "ListOfPairs" is undefined list of canonical pairs is NOT generated.
#undef ListOfPairs
//Equal to: Steps to form NonCanonical pairs are rejected and counted as iterations. (NCan_RC)

//If "BreakPairs" is defined, breaking pairs to form new pairs is allowed. (Paired_BC)
//#define BreakPairs

//If "BreakPairs" is undefined breaking pairs is not allowed. 
#undef BreakPairs
//Step that would potentially break a pair is rejected and counter as iteration. (Paired_RC)
//END STAS

using namespace std;

const double BoltzmannConstant = 0.001987; //constant in kcal/mol/K
const double delta_spring_const = 0;

struct NeigborDistance {
    int distance_left;
    int distance_right;
};


//This class provides an RNA class, and allows an array because there is a default constructor.
class RNA_Container {

public:

	//Default constructor
	RNA_Container();

	//Return a pointer to the start of the underlying RNA class.
	RNA *Return_RNA();

	//Allocate a member from the array of RNA classes.  Requires a null terminated c string with the sequence.
	//Also requires a specificationa as to whether this is RNA or DNA folding.
	//Return an int that indicates whether an error occurred while reading the file.  0 = no error, otherwise there was an error.
	int AllocateFromSequence(char *Sequence, bool IsRNA=true);


    //===> Allocate a member from the array of RNA classes. Requires the name of the structure.
    int AllocateFromStructure(char *Structure_Filename, bool IsRNA=true);


	//Destructor
	~RNA_Container();
	
private:

	//The underlying RNA class.
	RNA *rna;

};


//! MonteCarloClass Class.
/*!
	The MonteCarloClass class provides an entry point for using RNAstructure for Monte Carlo calculations.
*/

//Note the stylized comments provide facility for automatic documentation via doxygen.
class MonteCarloClass {

	//******************************
	//Constructors:
	//******************************
public:
	//!Constructor.

    //===> Initialize the MonteCarloClass with two structures
    MonteCarloClass(char *Structure1_Filename, char* Structure2_Filename, int replicas, double temp, double spring_constant, long seed=1, bool IsRNA=true);

	//!Constructor.

	//!	Provide a restart file name, a start temperature, a random seed, and whether the sequence is RNA or DNA.
	//!	This constructor generates internal error states that can be accessed by GetErrorState() after the constructor is called. 
	//!	\param Restart_Filename is a NULL terminated c string.  It is in the format written by WriteRestart.  
	//! \param seed is the temperature to start.
	//! \param seed is the random number seed for the random number generator.  Deafult is 1.
	//!	\param IsRNA is a bool that indicates whether this sequence is RNA or DNA.  true=RNA.  false=DNA.  Default is true.
	MonteCarloClass(char *Restart_Filename, double temp, double spring_constant, long seed=1, bool IsRNA=true);

	//!	Provide a restart file name, a random seed, and whether the sequence is RNA or DNA.
	//!	This constructor generates internal error states that can be accessed by GetErrorState() after the constructor is called. 
	//!	\param Restart_Filename is a NULL terminated c string.  It is in the format written by WriteRestart.  
	//! \param seed is the random number seed for the random number generator.  Deafult is 1.
	//!	\param IsRNA is a bool that indicates whether this sequence is RNA or DNA.  true=RNA.  false=DNA.  Default is true.
	MonteCarloClass(char *Restart_Filename, long seed=1, bool IsRNA=true);


	//! Build a list of all possible canonical pairs for move set.
	//! A pair list must be generated before starting the Monte Carlo steps.
	//! \return An int that gives an error state.  0 = no error, otherwise there was an error.  An error message can be retrieved using GetErrorMessage().
	int BuildPairList();

	//! Collect a monte carlo sample.

	//! Note that the ample is drawn on from replica zero, which is assumed to be the lowest temperature replica.
	//! \param interations is the number of iterations of sampling to perform.
	//! \param sampling_frequency is the frequency at which structure are sampled.
	//! \param ct_filename is a null terminated c string with the name of a ct file to which the sample is written.  This file is appended to, rather than overwritten.
	//! \return An int that gives an error state.  0 = no error, otherwise there was an error.  An error message can be retrieved using GetErrorMessage().
	// int sample(unsigned long iterations, unsigned int sampling_frequency, char *ct_filename);
	//===>Ju remove sampleling_frequency parameter
	int sample(unsigned long iterations, char *ct_filename);

	//! Set the name for the restart file.

	//! This filename is used for making a restart file.
	//! The constructor sets a default name of "MonteCarloRestart.ct".
	//! \param restartfilename is a pointer to cstring with the name of the file.
	void SetRestartFileName(char* restartfilename);

	void SetConstant(double spring_constant);

	void SetTemperatureRange(double low_temp, double high_temp);

	void SetWriteCTRange(unsigned long start_it, unsigned long end_it);

	//! Set the temperature of a replica.

	//! Indicate a replica and its temperature in Kelvins.
	//! Note that this temperature is for replica exhange, it does not affect the thermodynamics, which right now are fixed at 310.15 K.
	//! The constructor sets the temperatures of all replicas to 310.15 by default.
	//! Note that the code expects that the temperatures increase with increasing index.
	//! \param replica is an int that refers the replica number, which is zero-indexed.
	//! \param temperature is the temperature of a replica.
	//! \return An int that indiactes whether an error occurred.  0 = no error, otherwise an error occurred. An error message can be retrieved using GetErrorMessage(). 
	int SetReplicaTemperature(int replica, double temperature);

	
	//! Return an error message.

	//! If an error occurred, this function can return a meaningful error message.
	//! \return A string that explains the error that occurred.
	string GetErrorMessage();

	//! Return an error state, where a return of zero is no error.

	//!	This function returns an error flag.
	//!	An error of zero is no error.  Other returns indicate an error, where GetErrorMessage will provide a error message.
	//! \return An integer that provides the error state.
	int GetErrorState();
	
	//!Get the pointer to an underlying RNA class.

	//! Programmer provides an index to the replicas.  
	//! \param i is an int that is the replica being indexed.
	//! \return A pointer to an RNA class or NULL, if the int is out of range.
	RNA* GetRNAClass(int i);

	//! Write a restart file.

	//! Write a restart file that contains the current structures, and their temperatures.
	//! \param outputfilename is a pointer to cstring that specifies the name of the file to be written.
	//! \return An int that indicates the error state.  0 is no error.  (Note that no error handling has been yet programmed...).
	int WriteRestart(const char *outputfilename);
	
	//===>Ju
	//! This function calculate the distance.
	NeigborDistance Distance(RNA *rna1, RNA *rna2, RNA *rna3);
	//! This function calculate the tangent energy.
	double SpringEnergy(double spring_const_left, double spring_const_right, RNA *rna1, RNA *rna2, RNA *rna3);
	//===>Ju

	//!Destructor.
	~MonteCarloClass();

private:

	//A function to bundle the common elements of construction.
	void CommonConstructorElements();

	//Keep pointer to sequence filename.
	char *sequence_filename;



	//RNA_Container provides an array of RNA classes, to represent the replicas.
	RNA_Container *rna;

	//Keep track of the number of Replicas.
	int Replicas;

	//STAS: Uncomment this to generate list of canonical base pairs
#ifdef ListOfPairs
	//Keep track of possible moves. This is a list of pairs.  p5 has the 5' partner and p3 has the 3' partner, when the same index is used.
	vector<int> p5;
	vector<int> p3;
#endif
	//STAS: END

	//Set up new pair matrix:
	bool **pairs;

	//An array to keep track of the temperature of each replicate:
	double *temperature;

	//An array to keep track of the spring constant of the spring on the left of each image:
	double *springconst;

	//Keep track of energies:
	double *currentenergy;
	double *current_pseudo_energy;

	//Class to generate random numbers:
	randomnumber random;

	//Perform a step of replica exchange:
	void replicaexchange();

	//Perform one step of Monte Carlo:
	bool step();

	//This keeps track of whether an error occurred.
	int errorstate;

	//Keep track of what error occurred.
	string errorname;

	//The interval between replica exhange steps
	int replicaexchangerate;

	//The total number of replica exchanges
	int numreplicaexchanges;

    //The total number of accepted replica steps
    int TotalAcc;

    //The total number of rejected replica steps
    int TotalRej;


    // //The distantces of rna1 to rna2 (dist_a) and rna3 to rna2 (dist_b)
    // int dist_a, dist_b;
    
    //The tangent Energy.
    // double tangent_energy;
    // //===>Ju


// #ifdef REPLICASTATS
// 	//Holds the swap ststistics of replica exchange - which
// 	//replica swapped with which.
// 	vector<vector<int> > ReplicaStats;

// 	//Holds the swap ststistics of replica exchange - energy
// 	//difference of the swapped structures
// 	vector<vector<int> > ReplicaStatsEnergy;

//     //Holds the swap statistics of replica exchange - swap probability
//     vector<vector<double> > ReplicaStatsProb;

//     //Holds the swap statistics of replica exchange - random number
//     //of the swap
//     vector<vector<double> > ReplicaStatsTemproll;

//     //Holds the count of automatically accepted moves in MonteCarlo
//     vector<int> NucSwapStatsAutoAccepted;

//     //Holds the count of accepted moves in MonteCarlo after the random number roll
//     vector<int> NucSwapStatsRollAccepted;

//     //Holds the count of rejected moved in MonteCarlo after the random number roll
//     vector<int> NucSwapStatsRollRejected;

// #endif

	//The interval between writing restart files:
	int restartrate;

	//*notice
	int evaluaterate;
	//notice*

	//Keep a restart filename:
	string RestartFileName;

	//===>Ju
	//The spring constant
	double SpringConst; //constan in kcal/mol/distance

	//===>Ju
	//Setup temperature schedule points.
	double Temp;
	double LowTemp; 
	double HighTemp; 
	double DeltaTemp;

	//===>Ju
	//setup write ct range.
	unsigned long StartIt;
	unsigned long EndIt;

	//If this is replica exchange, allocate a pairlist array to facilite the swapping fo structures
	int *pairlist;

	//Vector to hold the order of replica exchanges
	vector<int> ReplicaOrder;


// #ifdef REPLICADEBUG
// 	//Name of file to output debuging data related to exchanges that are made.
// 	string debugingOutputName;

// 	//Count the number of exchanges occuring for development purposes
// 	int numberOfExchanges;
// 	ofstream outputFile;
// 	int *exchanges;
// 	int *accepted;
// 	int *rejected;
// 	int *pseudoknot;
// 	int *numIterations;
// 	int *alreadyOtherPairing;
// 	int *pairingMade;
// 	int *pairingBroken;
// 	vector<vector<int> > madeByPosition;
// 	vector<vector<int> > brokenByPosition;
// 	vector<vector<int> > chosenByPosition;
// 		#endif
};



#endif //!defined MONTECARLO_H
