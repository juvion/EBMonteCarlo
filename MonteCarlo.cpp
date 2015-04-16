// MonteCarlo.cpp : Defines the entry point for the console application//

#include "MonteCarloClass.h"


int main(int argc, char* argv[])
{
  	//An instance of the Monte Carlo class:
	MonteCarloClass *mc;
	//===>Ju: need set up a try error for iterations, should be compatible for temperature schedule.
	if (argc==13) {
		//This is an originating run from a sequence file:

		//Allocate the MonteCarloClass:
		//For now, use one replica, require one replica for now (1), require RNA folding (last true):
		mc = new MonteCarloClass(argv[1],argv[2],atoi(argv[6]),atof(argv[7]),atol(argv[5]),true);
		//check for errors:
		if (mc->GetErrorState()!=0) {
			//An error occurred:

			cerr << mc->GetErrorMessage().c_str() << "\n";
		}
		else {
			//No error yet:
			//Set the replica temperatures:
			//===>??? argv[7]>=2 ?
            // if(atoi(argv[7])!=1){
            if(atoi(argv[6])>=2){    
                //===> Set the temperature for the simulation
                mc->SetTemperatureRange(atof(argv[7]), atof(argv[8])); 
            }
            
			mc->SetRestartFileName(argv[9]);
			mc->SetConstant(atof(argv[10]));
			mc->SetWriteCTRange(atol(argv[11]), atol(argv[12]));

			//Build a list of allowed pairs:
			if (mc->BuildPairList()!=0) {
				//An error occurred
				cerr << mc->GetErrorMessage().c_str() << "\n";
			}
			else {
				//No error thus far, Actually do the sampling:
				if(mc->sample(atol(argv[4]),argv[3])!=0) {
					//An error occurred:
					cerr << mc->GetErrorMessage().c_str() << "\n";
				}
			}
		}

		delete mc;

	}
	else if (argc==11) {
		//This is a run from a restart file, but with temperature schedule:

		//Allocate the MonteCarloClass:
		//For now, use one replica, require one replica for now (1), require RNA folding (last true):
		mc = new MonteCarloClass(argv[1],atof(argv[5]),atol(argv[4]),true);
		//check for errors:
		if (mc->GetErrorState()!=0) {
			//An error occurred:

			cerr << mc->GetErrorMessage().c_str() << "\n";
		}
		else {
			//No error yet:
			//Set the replica temperatures:
            if(atoi(argv[6])>=2){    
                mc->SetTemperatureRange(atof(argv[5]), atof(argv[6])); 
            }
            
			mc->SetRestartFileName(argv[7]);
			mc->SetConstant(atof(argv[8]));
			mc->SetWriteCTRange(atol(argv[9]), atol(argv[10]));

			//Build a list of allowed pairs:
			if (mc->BuildPairList()!=0) {
				//An error occurred
				cerr << mc->GetErrorMessage().c_str() << "\n";
			}
			else {
				//No error thus far, Actually do the sampling:
				if(mc->sample(atol(argv[3]),argv[2])!=0) {
					//An error occurred:
					cerr << mc->GetErrorMessage().c_str() << "\n";
				}
			}
		}

		delete mc;
	}
	else if (argc==6) {
		//This is a restart run from a restart file (ct format):
		mc = new MonteCarloClass(argv[1],(long) atoi(argv[4]),true);
		//mc->SetReplicaExchangeFrequency(atoi(argv[6]));
		mc->SetRestartFileName(argv[5]);
		//check for errors:
		if (mc->GetErrorState()!=0) {
			//An error occurred:
			cerr << mc->GetErrorMessage().c_str() << "\n";
			
		}
		else {
			//No error on allocation:
			//Build a list of allowed pairs:
			if (mc->BuildPairList()!=0) {
				//An error occurred
				cerr << mc->GetErrorMessage().c_str() << "\n";

			}
			else {
				//No error thus far, Actually do the sampling:
				if(mc->sample(atol(argv[3]),argv[2])!=0) {
					//An error occurred:
					cerr << mc->GetErrorMessage().c_str() << "\n";

				}

			}

		}
		delete mc;

	}
	else {

		cout << "Usage: 0:[MonteCarlo] 1:[Input_structure1_name] 2:[Input_structure2_name] 3:[output_ct_filename] 4:[iterations] 5:[random_seed] 6:[number_of_replicas] 7:[lowest_replica_temperature(Kelvins)] 8:[highest_replica_temperature(Kelvins)] 9:[output_restartfilename] 10:[spring_constant] 11:[write_ct_start_it] 12:[write_ct_end_it]\n";
		cout << "or\n";
		cout << "Usage: 0:[MonteCarlo] 1:[input_restart_filename] 2:[output_ct_filename] 3:[iterations] 4:[random_seed] 5:[lowest_replica_temperature(Kelvins)] 6:[highest_replica_temperature(Kelvins)] 7:[output_restartfilename] 8:[spring_constant] 11:[write_ct_start_it] 9:[write_ct_end_it]\n";
		cout << "or\n";
		cout << "Usage: 0:[MonteCarlo] 1:[input_restart_filename] 2:[output_ct_filename] 3:[iterations] 4:[random_seed] 5:[output_restartfilename]\n";
				

	}

	return 0;
}

