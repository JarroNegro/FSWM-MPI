/**
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 */

#include <iostream> 
#include <stdlib.h>
#include <getopt.h>
#include <iomanip> 
#include "Sequence.h"
#include "Seed.h"

#include <mpi.h>
#include "tools.h"
 
void printHelp(){
 std::string help = 
    "\nUsage: ./fswm [options] <sequence> "
    "\n"
    "\n<sequence> format:"  
    "\n\t Sequence must be in FASTA format. All genomes must be contained in one FASTA file. Example:"
    "\n\t >Genome1"
    "\n\t ATAGTAGATGAT.."
    "\n\t >Genome2"
    "\n\t ATAGTAGATGAT.."
    "\n\t >Genome3"
    "\n\t ATGATGATGATGATG.."
    "\n\t .."
    "\n\t "
    "\nOptions:"        
    "\n\t -h: print this help and exit"
    "\n\t -k <integer>: pattern weight (default 12)"
    "\n\t -t <integer>: numer of threads (default: 10)"
    "\n\t -s <integer>: the minimum score of a spaced-word match to be considered homologous (default: 0)"
    "\n";
	std::cout << help << std::endl;
}

int weight = 12;
int dontCare = 100;
int threads = 10;
int threshold = 0;

std::string fileContent;

std::string patterns;




void parseParameters(int argc, char *argv[]){
	int option_char;
	 while ((option_char = getopt (argc, argv, "k:t:hs:")) != -1){ 
		switch (option_char){ 
			case 's': 
				threshold = atoi (optarg); 
				break;
			case 'k': 
				weight = atoi (optarg); 
				if(weight<8 || weight > 16){
					std::cerr << "Weight (-k) must be between 8 and 16"<< std::endl;
					exit (EXIT_FAILURE);
				}
				break;
			case 't': 
				threads = atoi (optarg); 
				if(threads<1){
					std::cerr << "threads (-t) must be an integer larger than 0"<< std::endl;
					exit (EXIT_FAILURE);
				}
				break;
			case 'h': 
				printHelp();
				exit (EXIT_SUCCESS);
				break;
			case '?': 
				printHelp();		
				exit (EXIT_FAILURE);
      	}
	}
}

void writeDmat(std::vector<std::vector<double> > dmat, std::vector<Sequence>& sequences){
	std::ofstream outfile;
	outfile.open("DMat");
	outfile << sequences.size() << std::endl;
	for (int i = 0; i < sequences.size(); i++) 
	{
		std::string name = sequences[i].getHeader();
		for(int k = 0; k < 10; k++){
			if(k >= name.length())
				outfile << " ";
			else
				outfile << name[k];
		}
     	for (int j = 0; j < sequences.size(); j++) 
     	{
			if (i > j) 
	    			outfile << std::fixed <<std::setprecision(12) << dmat[i][j] << "  ";
			else if(j>i)
				outfile << std::fixed<< std::setprecision(12) << dmat[j][i] << "  ";
			else
					outfile << std::setprecision(12) << "0" << "  ";
     	}
      		outfile << std::endl;
	}
	
	//std::cout<< dmat[1][0] << std::endl;
	outfile.close();
}

int main(int argc, char *argv[]){
	if(argc < 2)
	{
		printHelp();		
		exit (EXIT_FAILURE);
	}
	parseParameters(argc,  argv);
	std::string fileName(argv[argc-1]);

        MPI_Init(&argc, &argv);
        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);


	int n;
	std::vector<Sequence> sequences;


	int64_t len1;
        int len2;


	double inicio,fin,tiempo_reco_inicio,tiempo_reco_fin,tiempo_bcast_inicio,tiempo_bcast_fin;

	inicio=MPI_Wtime();


        Seed seed(weight,dontCare);
	std::vector<double> local_results;

	omp_set_dynamic(0);
        omp_set_num_threads(threads);


if(rank==0)
{

	//read file content
	fileContent= read_file(fileName);
        sequences = Sequence::read(fileContent);


	Seed::init();
        std::cout << sequences.size() << " sequences read Process:" <<rank << std::endl;

	if(sequences.size() < 2){
		std::cerr << "there must be at least 2 sequences"<< std::endl;
		exit (EXIT_FAILURE);
	}

        n = sequences.size();
        len1 = fileContent.size();

	patterns=join(seed.getPatterns());
        len2 = patterns.size();

//        std::cout << patterns <<  std::endl;

}

tiempo_bcast_inicio=MPI_Wtime();

        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&len1, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
        MPI_Bcast(&len2, 1, MPI_INT, 0, MPI_COMM_WORLD);


	if (rank != 0) {
	fileContent.resize(len1);
	patterns.resize(len2);
	}

        MPI_Bcast(&patterns[0], len2, MPI_CHAR, 0, MPI_COMM_WORLD);


    int64_t chunk_sizebcast = 2000000000;   // tamaño por bloque
    int64_t num_chunksbcast = (len1 + chunk_sizebcast - 1) / chunk_sizebcast;

for (int64_t chunkbcast = 0; chunkbcast < num_chunksbcast; ++chunkbcast) {
        int64_t startbcast = chunkbcast * chunk_sizebcast;
        int64_t endbcast = std::min(startbcast + chunk_sizebcast, len1);
        int current_chunk_sizebcast = endbcast - startbcast;



if(current_chunk_sizebcast>0)
{
       std::cout << " C:" << chunkbcast << " num:" << num_chunksbcast << " star:" << startbcast << " eNd:" << endbcast << " CURR:" << current_chunk_sizebcast << "\n"; 
       MPI_Bcast((void*)fileContent.c_str()+ startbcast , current_chunk_sizebcast, MPI_CHAR, 0, MPI_COMM_WORLD);

}


    }



tiempo_bcast_fin=MPI_Wtime();


        if (rank != 0) {

        Seed seed(weight,dontCare,split(patterns));
	seed.init();

         sequences.reserve(n);
         sequences = Sequence::read(fileContent);


 	// Seed::init();

	std::cout << sequences.size() << " sequences read Process:" <<rank << std::endl;
        if(sequences.size() < 2){
                std::cerr << "there must be at least 2 sequences"<< std::endl;
                exit (EXIT_FAILURE);
        }


        }


	int total_pairs = (n * (n - 1)) / 2;

	std::vector<std::pair<int, int>> pairs;
	for (int i = 0; i < n; ++i)
		for (int j = i + 1; j < n; ++j)
			pairs.emplace_back(i, j);

//	std::cout << "PAIRS:"<< pairs.size() << "\n";

	int chunk_size = total_pairs / size;
	int remainder = total_pairs % size;
	int start = rank * chunk_size + std::min(rank, remainder);
	int end = start + chunk_size + (rank < remainder ? 1 : 0);

//	std::cout << "Rank:" << rank << " Chunk:"<< chunk_size << " Remainder:"<< remainder << " Start:"<< start << " End:"<< end << "\n";
//	std::cout << "Proceso " << rank << " l:"<< seed.getLength() << " W:"<< seed.getWeight() << " DC:" << seed.getDontCare() << "\n";

	std::cout << "start sorting process:"<< rank << std::endl;
        #pragma omp parallel for schedule(runtime)
        for(int i = 0; i < sequences.size();i++)
        {
                if(sequences[i].getSequence().size()<1000){
                        std::cerr << "each sequence must be longer than 1000 base pairs"<< std::endl;
                        exit (EXIT_FAILURE);
                }
                sequences[i].sortFirstBits(seed);
                sequences[i].sortFirstBitsRev(seed);
        }

        for(int i = 0; i < sequences.size(); i++)
        {
                sequences[i].sortNextBits(seed);
                sequences[i].sortNextBitsRev(seed);
        }



   std::cout << "starting pairwise distance calculation "<<  sequences.size() << " Sequences Process:" << rank << std::endl;
	for (int idx = start; idx < end; ++idx) {
		int i = pairs[idx].first;
		int j = pairs[idx].second;

//	std::cout << "Calculating Process " << rank << " i:" << i << " j:" <<j  << "treads:" << threads << "Thres:" << threshold <<"\n";
	double result = sequences[i].compareSequences(sequences[j], seed, threads, threshold);
	local_results.push_back(result);
	}


tiempo_reco_inicio=MPI_Wtime();

if (rank != 0) {
	int count = local_results.size();
	MPI_Send(&count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	MPI_Send(local_results.data(), count, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
} else {
	std::vector<std::vector<double>> all_results(size);
	std::vector<std::vector<double> >DMat(sequences.size(), std::vector<double>(sequences.size(),0));
	all_results[0] = local_results;

	for (int src = 1; src < size; ++src) {
	int count;
	MPI_Recv(&count, 1, MPI_INT, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	all_results[src].resize(count);
	MPI_Recv(all_results[src].data(), count, MPI_DOUBLE, src, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	for (int p = 0; p < size; ++p) {

//		std::cout << "Resultados del proceso " << p << ":\n";
//std::cout << "Número de filas: " << DMat.size() << std::endl;
//if (!DMat.empty()) {
//std::cout << "Número de columnas (en la primera fila): " << DMat[0].size() << std::endl;
//}

		int start_p = p * chunk_size + std::min(p, remainder);
		int end_p = start_p + chunk_size + (p < remainder ? 1 : 0);
		int idx_p = start_p;
		for (double val : all_results[p]) {
		int i = pairs[idx_p].first;
                int j = pairs[idx_p].second;
//			std::cout << "i:" <<i << " j:"<<j << "V:" << val << " ";
                        DMat[i][j]=val;
                        DMat[j][i]=val;

		idx_p++;

		}
	}

	std::cout << std::endl << "done" << std::endl;
	writeDmat(DMat, sequences);
        std::cout << std::endl << "Distances written to file DMat" << std::endl;

}

tiempo_reco_fin=MPI_Wtime();

fin=MPI_Wtime();
std::cout << "Proceso:" <<rank << " Tiempo de broadcast:" << tiempo_bcast_fin - tiempo_bcast_inicio  << " Tiempo recolección:" << tiempo_reco_fin - tiempo_reco_inicio << " Tiempo total:"<< fin-inicio  << std::endl;

MPI_Finalize();
return 0;
}
