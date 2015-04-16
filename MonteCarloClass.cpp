#include "MonteCarloClass.h"
#include <iostream>
#include <sstream>
#include <math.h>
using namespace std;

//====>define constructor for RNA_Container Class.
RNA_Container::RNA_Container() {
    //set the pointer to RNA to NULL.
    rna = NULL;
}

//====>define Return_RNA function for RNA_Container Class.
RNA *RNA_Container::Return_RNA() {
    return rna;
}

//===> Initialize the RNA class using the sequence for RNA_Container Class.
int RNA_Container::AllocateFromSequence(char *Sequence, bool IsRNA) {
    rna = new RNA(Sequence,IsRNA);
    //Also specify a dummy pair to add the first structure
    rna->SpecifyPair(1,2);
    rna->RemoveBasePair(1);
    //Return any error codes from reading the file.
    return rna->GetErrorCode();
}

//===> Initialize the RNA class using the file name for the structure for RNA_Container Class.
int RNA_Container::AllocateFromStructure(char *Structure_Filename, bool IsRNA) {
    rna = new RNA(Structure_Filename,1,IsRNA);
    //Return any error codes from reading the file.
    return rna->GetErrorCode();
}

//===> define destructor for RNA_Container Class.
RNA_Container::~RNA_Container() {
    //If rna was allocated, delete it.
    if (rna!=NULL) {
        delete rna;
    }
}


//MonteCarloClass Constructor
MonteCarloClass::MonteCarloClass(char *Structure1_Filename, char* Structure2_Filename, int replicas, double temp, double spring_constant, long seed, bool IsRNA) {
    //Keep track of the number of Replicas.
    Replicas = replicas;

    //Determine the starting temperature:
    Temp = temp;

    //Allocate the temperature array:
    temperature = new double [Replicas];
    //Seed the random number generator:

    //notice*  
    //add a spring constant array for adaptive EB method.
    SpringConst = spring_constant;
    springconst = new double [Replicas - 1];
    //notice*//

    //===>random is a instance of random number class.
    random.seed(seed);

    //Keep track of the error state:
    errorstate = 0;

    //Keep an error message
    errorname = "No Errors.";

    //Allocate an array of RNA_Arrays.
    rna = new RNA_Container[Replicas];

    //For each element in rna, the RNA_Container array, allocate the underlying RNA class.
    //===> assign starting structure to the first half of Replicas
    for (int i=0;i<floor(Replicas/2);++i) {
        int error = rna[i].AllocateFromStructure(Structure1_Filename,IsRNA);
        if (error!=0) {
            //An error occurred:
            errorstate = 1;
            errorname = "Error in Allocate ";
            errorname +=  rna[i].Return_RNA()->GetErrorMessage(error);
        }
    }

    //===> assign ending structure to the other half of Replicas
    for (int i=floor(Replicas/2);i<Replicas;++i) {
        int error = rna[i].AllocateFromStructure(Structure2_Filename,IsRNA);
        if (error!=0) {
            //An error occurred:
            errorstate = 1;
            errorname = "Error in Allocate ";
            errorname +=  rna[i].Return_RNA()->GetErrorMessage(error);
        }
    }

    //Set the default temperature for all Replicas:
    for (int i=0; i<Replicas; ++i) {
        temperature[i] = Temp;
    }

    //*notice// set a spring constant array for adaptive EB method.
    //Set an array for spring constants:
    for (int i=0; i<Replicas-1; ++i) {
        springconst[i] = SpringConst;
    }
    //notice*//

    CommonConstructorElements();

    //create a vector of images to be shuffled.
    for(int i=1;i<Replicas-1;++i) {
        ReplicaOrder.push_back(i);
    }
}


//Alternative constructor from restart file using input structures only.
MonteCarloClass::MonteCarloClass(char *Restart_Filename, double temp, double spring_constant, long seed, bool IsRNA) {

    //Open the restart file in an insatnce of the RNA class:
    RNA *RestartFile;
    RestartFile = new RNA(Restart_Filename,1);

    //Determine the number of Replicas:
    Replicas = RestartFile->GetStructureNumber();

    //Determine the starting temperature:
    Temp = temp;

    //Allocate the temperature array:
    temperature = new double [Replicas];

    //notice*  
    //add a spring constant array for adaptive EB method.
    SpringConst = spring_constant;
    springconst = new double [Replicas - 1];
    //notice*//

    random.seed(seed);


    //Keep track of the error state:
    errorstate = 0;
    //Keep an error message
    errorname = "No Errors.";

    //Allocate an array of RNA_Arrays.
    rna = new RNA_Container[Replicas];

    //For each element in rna, the RNA_Container array, allocate the underlying RNA class.

    //allocate sequence to be large enough to hold the whole sequence and a null terminator:
    char *sequence;
    sequence = new char [RestartFile->GetSequenceLength()+1];

    //Now populate the sequence string with the nucleotides:
    for (int i=0;i<RestartFile->GetSequenceLength();++i) {
        sequence[i]=RestartFile->GetNucleotide(i+1);
    }
    sequence[RestartFile->GetSequenceLength()]='\0';

    for (int i=0;i<Replicas;++i) {
        //Set the default temperature for all Replicas:
        temperature[i] = Temp;

        int error = rna[i].AllocateFromSequence(sequence,IsRNA);
        if (error!=0) {
            //An error occurred:
            errorstate = 1;
            errorname = "Error in Allocation ";
            errorname +=  rna[i].Return_RNA()->GetErrorMessage(error);
        }
        else {
            //Now set the pairs:
            for (int j=1;j<=RestartFile->GetSequenceLength();++j) {

                if (RestartFile->GetPair(j,i+1)>j) {
                    rna[i].Return_RNA()->SpecifyPair(j,RestartFile->GetPair(j,i+1));
                }
            }
        }
    }

    //*notice// set a spring constant array for adaptive EB method.
    //Set an array for spring constants:
    for (int i=0; i<Replicas-1; ++i) {
        springconst[i] = SpringConst;
    }
    //notice*//


    //Clean up memory use:
    delete RestartFile;
    delete[] sequence;

    CommonConstructorElements();

    for(int i=1;i<Replicas-1;++i) {
        ReplicaOrder.push_back(i);
    }
}

//Alternative constructor from restart file. Restart with previous stopped temperature. It still needs to build previous temperature schedule.
MonteCarloClass::MonteCarloClass(char *Restart_Filename, long seed, bool IsRNA) {

    //Open the restart file in an instance of the RNA class:
    RNA *RestartFile;
    RestartFile = new RNA(Restart_Filename,1);

    //Determine the number of Replicas:
    Replicas = RestartFile->GetStructureNumber();

    //Allocate the temperature array:
    temperature = new double [Replicas];

    //*notice// add a spring constant array for adaptive EB method.
    //springconst = new double [Replicas - 1]
    //*notice//

    //Seed the random number generator:
    //===>random is a instance of random number class.
    random.seed(seed);


    //Keep track of the error state:
    errorstate = 0;
    //Keep an error message
    errorname = "No Errors.";

    //Allocate an array of RNA_Arrays.
    rna = new RNA_Container[Replicas];

    //For each element in rna, the RNA_Container array, allocate the underlying RNA class.

    //allocate sequence to be large enough to hold the whole sequence and a null terminator:
    char *sequence;
    sequence = new char [RestartFile->GetSequenceLength()+1];

    //Now populate the sequence string with the nucleotides:
    for (int i=0;i<RestartFile->GetSequenceLength();++i) {
        sequence[i]=RestartFile->GetNucleotide(i+1);
    }
    sequence[RestartFile->GetSequenceLength()]='\0';

    for (int i=0;i<Replicas;++i) {
        //Set the default temperature for all Replicas:
        string stringtemperature;
        stringtemperature = RestartFile->GetCommentString(i+1);
        cout << stringtemperature << endl;
        //stringtemperature.erase(stringtemperature.begin()+stringtemperature.find(' ',0),stringtemperature.end());
        temperature[i] = atof(stringtemperature.c_str());

        int error = rna[i].AllocateFromSequence(sequence,IsRNA);
        if (error!=0) {
            //An error occurred:
            errorstate = 1;
            errorname = "Error in Allocation ";
            errorname +=  rna[i].Return_RNA()->GetErrorMessage(error);
        }
        else {
            //Now set the pairs:
            for (int j=1;j<=RestartFile->GetSequenceLength();++j) {

                if (RestartFile->GetPair(j,i+1)>j) {
                    rna[i].Return_RNA()->SpecifyPair(j,RestartFile->GetPair(j,i+1));
                }
            }
        }
    }

    //*notice// set a spring constant array for adaptive EB method.
    //Set an array for spring constants, read in spring constant from restart:
    // for (int i=0; i<Replicas-1; ++i) {
    //     springconst[i] = SpringConst;
    // }
    //*notice//

    //Clean up memory use:
    delete RestartFile;
    delete[] sequence;

    CommonConstructorElements();

    for(int i=1;i<Replicas-1;++i) {
        ReplicaOrder.push_back(i);
    }
}

//canonical base pair rule.
inline bool allowedpair(char n1, char n2) {
    if (n1=='A'&&n2=='U') return true;
    else if (n1=='C'&&n2=='G') return true;
    else if (n1=='G'&&(n2=='C'||n2=='U')) return true;
    else if (n1=='U'&&(n2=='A'||n2=='G')) return true;
    else return false;
}

// Build a list of all possible canonical pairs for move set.
//Currently allow any pair.
int MonteCarloClass::BuildPairList() {
    int i,j;

    pairs = new bool *[rna[0].Return_RNA()->GetSequenceLength()+1];
    for (i=1;i<=rna[0].Return_RNA()->GetSequenceLength();++i) {
        pairs[i] = new bool [rna[0].Return_RNA()->GetSequenceLength()+1];

        for (j=1;j<=rna[0].Return_RNA()->GetSequenceLength();++j) {
            pairs[i][j]=false;
        }
    }

    //2-D search across all possible pairing partners:
    for (i=1;i<=rna[0].Return_RNA()->GetSequenceLength();++i) {
        for (j=i+minloop+1;j<=rna[0].Return_RNA()->GetSequenceLength();++j) {

            if (allowedpair(rna[0].Return_RNA()->GetNucleotide(i),rna[0].Return_RNA()->GetNucleotide(j))) {
                //also, try making sure this is not an isolated base pair
                bool before = false;
                if ((i>1&&j!=rna[0].Return_RNA()->GetSequenceLength())) {
                    before = allowedpair(rna[0].Return_RNA()->GetNucleotide(i-1),rna[0].Return_RNA()->GetNucleotide(j+1));
                    //before = inc[ct->numseq[i-1]][ct->numseq[j+1]];
                }

                bool after;
                //after = 0 if a stacked pair cannot form 3' to i
                if (((j-i)>minloop+2)&&(i!=rna[0].Return_RNA()->GetSequenceLength())) {
                    after = allowedpair(rna[0].Return_RNA()->GetNucleotide(i+1),rna[0].Return_RNA()->GetNucleotide(j-1));//inc[ct->numseq[i+1]][ct->numseq[j-1]];
                }

                else after = false;

                //if there is a stackable pairs to i.j then add that pair to the list
                if ((before)||(after)) {
                    //STAS: Uncomment this code to generate list of canonical base pairs
#ifdef ListOfPairs
                    p5.push_back(i);
                    p3.push_back(j);
#endif
                    //END STAS

                    pairs[i][j]=true;
                }
            }
        }
    }

    //Right now, there is no error checking, so return that no error occurred:
    return 0;
}



//Perform the Monte Carlo sampling.
int MonteCarloClass::sample(unsigned long iterations, char *ct_filename) {

    //calculate for first and last image's free energy.
    currentenergy[0] = rna[0].Return_RNA()->CalculateFreeEnergy(1,false);
    currentenergy[Replicas-1] = rna[Replicas-1].Return_RNA()->CalculateFreeEnergy(1,false);

    current_pseudo_energy[0] = 0;
    current_pseudo_energy[Replicas-1] = 0;

    //calculate for intermediate images' structures free energies.
    for (int rep=1; rep<Replicas-1; ++rep) {
        currentenergy[rep] = rna[rep].Return_RNA()->CalculateFreeEnergy(1,false);
        int errorCode1 = rna[rep].Return_RNA()->GetErrorCode();
        if(errorCode1 != 0) {
                //An error occurred, so store a string with a message an return 1:
                errorname = "Error in CalculateFreeEnergy:";
                errorname+= rna[rep].Return_RNA()->GetErrorMessageString(errorCode1);
                return 1;
        }
    }

    //initialize current pseudo energies
    for (int rep=1; rep<Replicas-1; ++rep) {
        current_pseudo_energy[rep] = SpringEnergy(springconst[rep-1], springconst[rep], rna[rep - 1].Return_RNA(), rna[rep].Return_RNA(), rna[rep + 1].Return_RNA()) + 
                                     sqrt((currentenergy[rep] - currentenergy[rep-1]) * (currentenergy[rep] - currentenergy[rep-1]) / 2 + 
                                         (currentenergy[rep+1] - currentenergy[rep]) * (currentenergy[rep+1] - currentenergy[rep]) / 2);
        int errorCode1 = rna[rep].Return_RNA()->GetErrorCode();
        if(errorCode1 != 0) {
                //An error occurred, so store a string with a message an return 1:
                errorname = "Error in CalculateFreeEnergy:";
                errorname+= rna[rep].Return_RNA()->GetErrorMessageString(errorCode1);
                return 1;
        }
    }


    //Now perform the sampling
    //setup temperature schedule
    int steps_stage1 = iterations / 5;
    int steps_stage2 = iterations / 4;
    int steps_stage3 = iterations - steps_stage1 - steps_stage2;
    double increase_delta = (HighTemp - LowTemp) / double(steps_stage1);
    double decrease_delta = (HighTemp - LowTemp) / double(steps_stage3);

    //start the iteration of simulation
    for (unsigned long it=0; it < iterations; ++it) {

        //STAS: Uncomment the following line to enable the counter
#ifdef ItCounter
        if((it+1) % 100000 == 0) cout << "\r" << it+1 << flush;
#endif
        //END STAS
        

        //*notice //evaluate the distance between two neighbors and adjust the spring constant for every evaluaterate steps.
        if ((it%evaluaterate==0)&&(it> 0)) {
            for (int rep=1; rep<Replicas-1; ++rep) {
                NeigborDistance neighbor_distances = Distance(rna[rep - 1].Return_RNA(), rna[rep].Return_RNA(), rna[rep + 1].Return_RNA());
                int dist_left = neighbor_distances.distance_left;
                // int dist_right = neighbor_distances.distance_right;

                if (dist_left > 1) {
                    springconst[rep-1] += delta_spring_const; //delta_spring_const = 0.01
                } else if (springconst[rep-1] - delta_spring_const > SpringConst && dist_left < 1) {
                    springconst[rep-1] -= delta_spring_const; 
                }
                // if (dist_right > 1){
                //     springconst[rep] += delta_spring_const;
                // } else if (springconst[rep] - delta_spring_const > SpringConst) {
                //     springconst[rep] -= delta_spring_const;
                // }
            } 
        }
        //notice*  

        //Write a restart file, if needed:
        if ((it%restartrate==0)&&(it>0)) {
            WriteRestart(RestartFileName.c_str());
        }

        if (it < steps_stage1) {
            for (int rep=0; rep<Replicas; ++rep) temperature[rep] = LowTemp;
            LowTemp += increase_delta;
        }

        else if (it >= steps_stage1 && it < (steps_stage1 + steps_stage2)) {
            for (int rep=0; rep<Replicas; ++rep) temperature[rep] = HighTemp;
        }

        else {
            HighTemp -= decrease_delta;
            for (int rep=0; rep<Replicas; ++rep) temperature[rep] = HighTemp;
        }
        //execute the step.
        bool is_step_accept = step();
    }

    //Always write a restart file at the end of the run:
    WriteRestart(RestartFileName.c_str());

    ///STAS changed it to write new .ct file insteas of appending
    rna[0].Return_RNA()->WriteCt(ct_filename, false);
    for (int i = 1; i < Replicas; i++) {
        rna[i].Return_RNA()->CalculateFreeEnergy(1,false);
        rna[i].Return_RNA()->WriteCt(ct_filename, true);
    }
    //No errors found, return 0:
    return 0;
}


//Set the name of the restart file:
void MonteCarloClass::SetRestartFileName(char* restartfilename) {

    RestartFileName = restartfilename;

}

void MonteCarloClass::SetConstant(double spring_constant) {

    SpringConst = spring_constant;

}


void MonteCarloClass::SetTemperatureRange(double low_temp, double high_temp) {

    LowTemp = low_temp;
    HighTemp = high_temp;

}

void MonteCarloClass::SetWriteCTRange(unsigned long start_it, unsigned long end_it) {

    StartIt = start_it;
    EndIt = end_it;

}

//Set the temperature of a replica:
int MonteCarloClass::SetReplicaTemperature(int replica, double replica_temperature) {

    if (replica<0||replica>=Replicas) {
        //This index does not make sense.
        errorname = "Index out of range.";
        return 1;
    }

    //Past error check:
    temperature[replica] = replica_temperature;
    return 0;

}

//Return an error message as a string.
string MonteCarloClass::GetErrorMessage() {

    return errorname;

}

// Return an error state, where a return of zero is no error.
int MonteCarloClass::GetErrorState() {

    return errorstate;

}

//Write a restart file.
int MonteCarloClass::WriteRestart(const char *outputfilename) {

    char *sequence;

    //allocate sequence to be large enough to hold the whole sequence and a null terminator:
    sequence = new char [rna[0].Return_RNA()->GetSequenceLength()+1];

    //Now populate the sequence string with the nucleotides:
    for (int i=0;i<rna[0].Return_RNA()->GetSequenceLength();++i) {

        sequence[i]=rna[0].Return_RNA()->GetNucleotide(i+1);

    }
    sequence[rna[0].Return_RNA()->GetSequenceLength()]='\0';

    //Use an RNA class to contain all the structures:
    RNA *restartcollector;
    restartcollector = new RNA(sequence);

    //Now populate the structures:
    for (int i =0; i< Replicas; ++i) {

        //Specify a dummy pair.  This ensures the underlying memory is correctly allocated for the structures.
        restartcollector->SpecifyPair(1,2,i+1);
        restartcollector->RemoveBasePair(1,i+1);

        for (int j = 1 ; j <= rna[0].Return_RNA()->GetSequenceLength(); ++j) {

            if (rna[i].Return_RNA()->GetPair(j)>j) {
                restartcollector->SpecifyPair(j,rna[i].Return_RNA()->GetPair(j),i+1);
            }
        }

        //Put the temperature at the start of the structure comment
        string comment;
        char commentchar[20];
        sprintf(commentchar,"%f",temperature[i]);
        comment = commentchar;
        comment += " ";
        //comment += rna[0].Return_RNA()->GetCommentString();
        restartcollector->AddComment(comment.c_str(),i+1);
    }

    //*notice//
    //write a loop to record spring constant in the restart file head, with temperature.
    //*notice//


    restartcollector->WriteCt(outputfilename);

    //cleanup memory use:
    delete[] sequence;
    delete restartcollector;

    return 0;//For now, there is no error trapping.
}

//Destructor
MonteCarloClass::~MonteCarloClass() {

    if (pairs!=NULL) {

        for (int i=1;i<=rna[0].Return_RNA()->GetSequenceLength();++i) delete[] pairs[i];
        delete[] pairs;

    }

    //Clean up memory use.

    //These were allocated in the Constructor:
    delete[] rna;
    delete[] temperature;
    delete[] springconst;
    delete[] currentenergy;
    delete[] current_pseudo_energy;
 

    //If this is replica exchange, clean up the pairlist array
    if (Replicas>1) {

        delete[] pairlist;

    }
}

//*notice//
//spring constant array may be created inside this function
//*notice//
void MonteCarloClass::CommonConstructorElements() {

    //Keep track of structure folding free energies
    currentenergy = new double [Replicas];
    current_pseudo_energy = new double [Replicas];

    //If this is replica exchange, allocate the pairlist array to facilitate the swapping of pair information

    if (Replicas>1) {
        pairlist = new int [rna[0].Return_RNA()->GetSequenceLength()+1];
    }

    //Set a default rate for writing to the restart file:
    restartrate = 10000;

    //*notice
    evaluaterate = 100; 
    //notice*

    //Set a default filename for writing restart files:
    RestartFileName = "MonteCarloRestart.ct";

    pairs=NULL;

}



//*notice
//calculate the distances between one image's two neighbors.
NeigborDistance MonteCarloClass::Distance(RNA *rna1, RNA *rna2, RNA *rna3) {
    // NeigborDistance dist;
    NeigborDistance neighbor_distances;
    int dist_a = 0;
    int dist_b = 0;
    //iterate the each bases of the sequence.
    for (int i = 1; i <= rna1->GetSequenceLength() - 1; i++) {
        for (int j = i + 1; j <= rna1->GetSequenceLength(); j++) {
            //find the pair to base index for rna1, rna2 and rna3
            int pair_to1 = rna1->GetPair(i);
            int pair_to2 = rna2->GetPair(i);
            int pair_to3 = rna3->GetPair(i);
            //if pair(i, j) exists in rna2.
            if (pair_to2 == j) {
                //if pair(i, j) does not exist in rna1, dist_a increment 2.
                if (pair_to1 != j) dist_a += 1;
                //if pair(i, j) does not exist in rna3, dist_b increment 2.
                if (pair_to3 != j) dist_b += 1;
            }

            //if pair(i, j) does not exist in rna2.
            else {
                //if pair(i, j) exists in rna1, dist_a increment 2.
                if (pair_to1 == j) dist_a += 1;
                //if pair(i, j) exists in rna3, dist_b increment 2.
                if (pair_to3 == j)  dist_b += 1;
            }
        }
    }
    //calculate the distance
    // std::cout << dist_a << ' ' << dist_b << std::endl;
    neighbor_distances.distance_left = dist_a;
    neighbor_distances.distance_right = dist_b;
    return neighbor_distances;
}
//notice*

//*notice//
//other methods for spring energy calculation could be found from previous Monte Carlo code. 
//notice*//
//tangentEnergy2.6.1, counting basepairs not bases. The sprint constant needs to be 4 times larger
// linear search method to search base, and the result is same with tangentEnergy1.6
double MonteCarloClass::SpringEnergy(double spring_const_left, double spring_const_right, RNA *rna1, RNA *rna2, RNA *rna3) {
    int dist_a = 0;
    int dist_b = 0;
    //iterate the each bases of the sequence.
    for (int i = 1; i <= rna1->GetSequenceLength() - 1; i++) {
        for (int j = i + 1; j <= rna1->GetSequenceLength(); j++) {
            //find the pair to base index for rna1, rna2 and rna3
            int pair_to1 = rna1->GetPair(i);
            int pair_to2 = rna2->GetPair(i);
            int pair_to3 = rna3->GetPair(i);
            //if pair(i, j) exists in rna2.
            if (pair_to2 == j) {
                //if pair(i, j) does not exist in rna1, dist_a increment 2.
                if (pair_to1 != j) dist_a += 1;
                //if pair(i, j) does not exist in rna3, dist_b increment 2.
                if (pair_to3 != j) dist_b += 1;
            }

            //if pair(i, j) does not exist in rna2.
            else {
                //if pair(i, j) exists in rna1, dist_a increment 2.
                if (pair_to1 == j) dist_a += 1;
                //if pair(i, j) exists in rna3, dist_b increment 2.
                if (pair_to3 == j)  dist_b += 1;
            }
        }
    }

    //*notice//
    //with array of spring constant, the equation to calculate the energy changes. 
    //calculate the spring_energy.
    // double spring_energy = SpringConst * (double(dist_a * dist_a) + double(dist_b * dist_b));
    double spring_energy = spring_const_left * (double(dist_a * dist_a)) + spring_const_right * (double(dist_b * dist_b));
    return spring_energy;
    //notice*//
}


//Perform one step of Monte Carlo:
bool MonteCarloClass::step() {
    int *rollon5, *rollon3;
    double *randomarray;
    //by default, is_step_accept is set to false. If any image's new move is accepted,
    //then it is changed to true.
    bool is_step_accept = false;

    // std::cout << proceed << std::endl;
    //Shuffle the ReplicaOrder [Replicas] times
    for(int i=0;i<ReplicaOrder.size();++i){
        //Temporarily store the current "order" at position i
        int tempOrder = ReplicaOrder[i];
        //Run the random number generator to get the position for swapping and store it in 'swap'
        //===>Ju changed the shuffle process exclude the first and last replica.
        // int swapTo = (int)floor((double)(Replicas-1)*random.roll());
        int swapTo = (int)floor((double)(ReplicaOrder.size())*random.roll());
        //Swap the "order" from position determined by the random number to position 'i'
        ReplicaOrder[i]=ReplicaOrder[swapTo];
        //Assign the original "order" from position 'i' to the position determined by the random number
        ReplicaOrder[swapTo]=tempOrder;
    }  

    rollon5 = new int[Replicas];
    rollon3 = new int[Replicas];
    randomarray = new double[Replicas];

    //====> image[1] to image[Replicas-2] need to be calculated.
    // for (int rep=0; rep<Replicas; ++rep) {
    for (int rep=1; rep<Replicas-1; ++rep) {
        //randomly choose a 5' nuc and 3' nuc
        rollon5[rep]=1+(int) (random.roll()*((double) (rna[rep].Return_RNA()->GetSequenceLength())));
        rollon3[rep]=1+(int) (random.roll()*((double) (rna[rep].Return_RNA()->GetSequenceLength())));
        //randomly choose a number to compare to Boltzmann prob.
        randomarray[rep]=random.roll();
    }

    //Do this for each image, the order was shuffled:
    for (int i=0; i<ReplicaOrder.size(); ++i) {
        bool proceed = true;
        int rep = ReplicaOrder[i];

        //*notice //this energy calculation may be redundant or may need to move to other place to increase the efficiency.
        //calculate the current spring energy and free energy before the new base pair is created or disrupted.
        // double current_free_energy = rna[rep].Return_RNA()->CalculateFreeEnergy(1,false);
        // double current_spring_energy = SpringEnergy(springconst[rep], springconst[rep], rna[rep - 1].Return_RNA(), rna[rep].Return_RNA(), rna[rep + 1].Return_RNA());
        //notice*//

        //Generate the nucleotides using random numbers.
        int on5 = rollon5[rep];
        int on3 = rollon3[rep];
        //exchange the 5' and 3' nucleotides if 3'-nt is at 5' side. index of nt is from 5' to 3'
        if (on3<on5) {
            int tmp = on5;
            on5=on3;
            on3=tmp;
        }

        //Store the original values
        int originalnuc1 = rna[rep].Return_RNA()->GetPair(on5);
        int error5 = rna[rep].Return_RNA()->GetErrorCode();
        if(error5 !=0) {
            cout<< "Error in GetPair:" <<rna[rep].Return_RNA()->GetErrorMessageString(error5);
        }
        int originalnuc2 = rna[rep].Return_RNA()->GetPair(on3);
        int error6 = rna[rep].Return_RNA()->GetErrorCode();
        if(error6 !=0) {
            cout<< "Error in GetPair:"<<rna[rep].Return_RNA()->GetErrorMessageString(error6);
        }

        //case1: Now see if the nucs are paired to each other, if so, break the pair
        if (originalnuc1==on3) {
            //break the pair of nuc1 and nuc2
            int error3 = rna[rep].Return_RNA()->SpecifyPair(originalnuc1,0);
            int error4 = rna[rep].Return_RNA()->SpecifyPair(originalnuc2,0);
            // std::cout << on3 << ' ' << on5 << ' ' << proceed << std::endl;

            if(error3!=0) {
                cout << "Error in SpecifyPair:"<<rna[rep].Return_RNA()->GetErrorMessageString(error3);
            }
            if(error4!=0) {
                cout << "Error in SpecifyPair:"<<rna[rep].Return_RNA()->GetErrorMessageString(error4);
            }
        }

        //case2: if originalnuc1 and originalnuc2 paired to other bases, not to each other, then NOT proceed.
        else if (originalnuc1!=0 || originalnuc2!=0) {
                proceed = false;
        }
        //case3: on5 and on3 are not originally paired to any nucleotide.
        else {
            //only when on5 and on3 are able to form canonical pair, then proceed.
            proceed = pairs[on5][on3];
            //if on5 and on3 can form canonical pair, check for pseudoknot.
            if (proceed){
                //Now check for pseudoknot
                for (int j=1; j<on5 && proceed; j++) {
                    //Found pk, then not proceed
                    if (rna[rep].Return_RNA()->GetPair(j)>on5&&rna[rep].Return_RNA()->GetPair(j)<on3){
                        proceed=false;
                    }
                }
                for (int j=on5+1; j<on3 && proceed; j++) {
                    //Found pk, then not proceed
                    if (rna[rep].Return_RNA()->GetPair(j)>on3){
                        proceed=false;
                    }
                }
            }
            if (proceed) {
                //No pseudoknot found, go ahead and specify the pairs:
                // rna[rep].Return_RNA()->SpecifyPair(originalnuc1,0);//Use only if breaking pairs to form new pairs is allowed.
                // rna[rep].Return_RNA()->SpecifyPair(originalnuc2,0);//Use only if breaking pairs to form new pairs is allowed.
                rna[rep].Return_RNA()->SpecifyPair(on5,on3);
            }
        }


//################# passed checks: 1. base pairing or breaking only upon canonical pairs; 2. no pseudoknot. #############
//################# process new checks:  #############
//         if (proceed) {
//             //calculate free energy, spring energy and total energy for new conformation.
//             double proceed_spring_energy = SpringEnergy(springconst[rep-1], springconst[rep], rna[rep - 1].Return_RNA(), rna[rep].Return_RNA(), rna[rep + 1].Return_RNA());
//             double proceed_energy = rna[rep].Return_RNA()->CalculateFreeEnergy(1,false) + proceed_spring_energy;// changed to false, since new energy model is used.
//             //if proceed_energy is higher than current energy, then perform a Boltzmann test. 
//             if (proceed_energy > currentenergy[rep]) {
//                 //########################################
//                 //# ====> calculate probability based on Boltzmann Distribution. Equation
//                 //# F(state) is proportional to exp(-E/kT), E is the energy, k is Boltzmann Constant, T is temperature.
//                 //# F(state2)/F(state1) = exp((E1 - E2)/kT)
//                 //#########################################
//                 double prob = std::exp((currentenergy[rep] - proceed_energy)/(BoltzmannConstant*temperature[rep]));
//                 if (randomarray[rep]>prob) {
//                     //Don't accept the structure, revert back to the original
//                     rna[rep].Return_RNA()->SpecifyPair(originalnuc1,on5);
//                     rna[rep].Return_RNA()->SpecifyPair(originalnuc2,on3);
//                     //update the folding free energy, since it was recalculated previously with new proposed structure.
//                     rna[rep].Return_RNA()->CalculateFreeEnergy(1,false);
//                     proceed = false;
//                 }
//                 else {
//                     // cout << temperature[rep] << endl;
//                     //accept the structure, store the energy, and update the neighbor images' current energy, since spring energy changed.
//                     currentenergy[rep] = proceed_energy;
//                     //for image1 and image0, no left neighbor energy is needed to update
//                     if (rep > 1) currentenergy[rep-1] = rna[rep-1].Return_RNA()->GetFreeEnergy(1) + SpringEnergy(springconst[rep-2], springconst[rep-1], rna[rep - 2].Return_RNA(), rna[rep-1].Return_RNA(), rna[rep].Return_RNA());
//                     //for image[Replicas-1] and image[Replicas-2] no right neighbor energy is needed to update
//                     if (rep < Replicas-2 ) currentenergy[rep+1] = rna[rep+1].Return_RNA()->GetFreeEnergy(1) + SpringEnergy(springconst[rep], springconst[rep+1], rna[rep].Return_RNA(), rna[rep+1].Return_RNA(), rna[rep+2].Return_RNA());
//                 }
//             }

//             //if proceed_energy is lower than or equal to current energy, then accept the new conformation.
//             else {
//                 //accept the structure, so store the energy, and update the neighbor images' current energy, since spring energy changed.
//                 currentenergy[rep] = proceed_energy;
//                 if (rep > 1) currentenergy[rep-1] = rna[rep-1].Return_RNA()->GetFreeEnergy(1) + SpringEnergy(springconst[rep-2], springconst[rep-1], rna[rep - 2].Return_RNA(), rna[rep-1].Return_RNA(), rna[rep].Return_RNA());
//                 //for image[Replicas-1] and image[Replicas-2] no right neighbor energy is needed to update
//                 if (rep < Replicas-2 ) currentenergy[rep+1] = rna[rep+1].Return_RNA()->GetFreeEnergy(1) + SpringEnergy(springconst[rep], springconst[rep+1], rna[rep].Return_RNA(), rna[rep+1].Return_RNA(), rna[rep+2].Return_RNA());
//             }
//         }
//     is_step_accept = is_step_accept || proceed;
//     }//end of a rep

//     delete[] rollon5;
//     delete[] rollon3;
//     delete[] randomarray;
//     return is_step_accept;
// }
       if (proceed) {
            //calculate root mean of the square of free energy gaps between neighbors. sseg = (E[rep] - E[rep-1])^2 + (E[rep+1] - E[rep])^2
            double proceed_free_energy = rna[rep].Return_RNA()->CalculateFreeEnergy(1,false);
            double proceed_spring_energy = SpringEnergy(springconst[rep-1], springconst[rep], 
                                                        rna[rep - 1].Return_RNA(), rna[rep].Return_RNA(), rna[rep + 1].Return_RNA());

            // double proceed_rms_dG = sqrt((proceed_free_energy - currentenergy[rep-1]) * (proceed_free_energy - currentenergy[rep-1]) / 2 + 
            //                        (currentenergy[rep+1] - proceed_free_energy) * (currentenergy[rep+1] - proceed_free_energy) / 2)

            //compute the pseudo energy as sum of proceed_spring_energy and proceed_rms_dG
            double proceed_pseudo_energy = proceed_spring_energy + 
                                           sqrt((proceed_free_energy - currentenergy[rep-1]) * (proceed_free_energy - currentenergy[rep-1]) / 2 + 
                                           (currentenergy[rep+1] - proceed_free_energy) * (currentenergy[rep+1] - proceed_free_energy) / 2);

            //if proceed pseudo energy is lower than previous conformation, then accept the conformation.
            if (proceed_pseudo_energy < current_pseudo_energy[rep]) {
                //accept the structure, and store the energy, and update the neighbor images' current energy, since spring energy changed.
                currentenergy[rep] = proceed_free_energy;
                current_pseudo_energy[rep] = proceed_pseudo_energy;
                if (rep > 1) current_pseudo_energy[rep-1] = sqrt((currentenergy[rep-1] - currentenergy[rep-2]) * (currentenergy[rep-1] - currentenergy[rep-2]) / 2 + 
                                                            (currentenergy[rep] - currentenergy[rep-1]) * (currentenergy[rep] - currentenergy[rep-1]) / 2) + 
                                                            SpringEnergy(springconst[rep-2], springconst[rep-1], rna[rep - 2].Return_RNA(), rna[rep-1].Return_RNA(), rna[rep].Return_RNA());
                //for image[Replicas-1] and image[Replicas-2] no right neighbor energy is needed to update
                if (rep < Replicas-2 ) current_pseudo_energy[rep+1] = sqrt((currentenergy[rep+1] - currentenergy[rep]) * (currentenergy[rep+1] - currentenergy[rep]) / 2 + 
                                                            (currentenergy[rep+2] - currentenergy[rep+1]) * (currentenergy[rep+2] - currentenergy[rep+1]) / 2) + 
                                                            SpringEnergy(springconst[rep], springconst[rep+1], rna[rep].Return_RNA(), rna[rep+1].Return_RNA(), rna[rep+2].Return_RNA());
            }
            else {
                //################Apply Boltzmann check on pseudo energy#######
                //# ====> calculate probability based on Boltzmann Distribution. Equation
                //# F(state) is proportional to exp(-E/kT), E is the energy, k is Boltzmann Constant, T is temperature.
                //# F(state2)/F(state1) = exp((E1 - E2)/kT)
                //#############################################################
                double prob = std::exp((current_pseudo_energy[rep] - proceed_pseudo_energy)/(BoltzmannConstant*temperature[rep])); 

                if (randomarray[rep]>prob) {
                    //Don't accept the structure, revert back to the original
                    rna[rep].Return_RNA()->SpecifyPair(originalnuc1,on5);
                    rna[rep].Return_RNA()->SpecifyPair(originalnuc2,on3);
                    //update the folding free energy, since it was recalculated previously with new proposed structure.
                    rna[rep].Return_RNA()->CalculateFreeEnergy(1,false);
                    proceed = false;
                }
                else{
                    //accept the structure, and store the energy, and update the neighbor images' current energy, since spring energy changed.
                    currentenergy[rep] = proceed_free_energy;
                    current_pseudo_energy[rep] = proceed_pseudo_energy;
                    if (rep > 1) current_pseudo_energy[rep-1] = sqrt((currentenergy[rep-1] - currentenergy[rep-2]) * (currentenergy[rep-1] - currentenergy[rep-2]) / 2 + 
                                                                (currentenergy[rep] - currentenergy[rep-1]) * (currentenergy[rep] - currentenergy[rep-1]) / 2) + 
                                                                SpringEnergy(springconst[rep-2], springconst[rep-1], rna[rep - 2].Return_RNA(), rna[rep-1].Return_RNA(), rna[rep].Return_RNA());
                    //for image[Replicas-1] and image[Replicas-2] no right neighbor energy is needed to update
                    if (rep < Replicas-2 ) current_pseudo_energy[rep+1] = sqrt((currentenergy[rep+1] - currentenergy[rep]) * (currentenergy[rep+1] - currentenergy[rep]) / 2 + 
                                                                (currentenergy[rep+2] - currentenergy[rep+1]) * (currentenergy[rep+2] - currentenergy[rep+1]) / 2) + 
                                                                SpringEnergy(springconst[rep], springconst[rep+1], rna[rep].Return_RNA(), rna[rep+1].Return_RNA(), rna[rep+2].Return_RNA());                    
                }

            }
        }
    is_step_accept = is_step_accept || proceed;
    }//end of a rep

    delete[] rollon5;
    delete[] rollon3;
    delete[] randomarray;
    return is_step_accept;
}