#include "Algorithm.h"
#define max(x,y) ( x < y ? y : x )

Algorithm::Algorithm(int R, int Q, ofstream * pfichier){
	
	penalExtension = R;
	penalOpen = Q;

	// écrire les informations de l'algorithme dans le fichier texte
	if (pfichier->is_open()){
		*pfichier << "ALGORITHM INFORMATION :" << endl;
		*pfichier << "Score matrix = BLOSSUM" << endl;
		*pfichier << "Extension Penalty = " << R << endl;
		*pfichier << "Open Penalty = " << Q << endl << endl;
	}
}

int Algorithm::calculateScore(int * querySequence, int sequenceSize, int sequenceSize2, int m_matrixBlosum[28][28], int * databaseSequence){
	
	// calculer le score d'un alignement entre deux séquences, basé sur la matrice blosum
	
	int* matrixH = new int[sequenceSize+1];
	int* matrixE = new int[sequenceSize+1];
	int F = 0;
	int H_up = 0;
	int H_diag = 0;
	int H = 0;
	
	int scoreMax =0;
	
	for (unsigned int i=0; i< sequenceSize+1;i++)
		{
		matrixE[i] = 0;
		matrixH[i] = 0;
		}

	for (unsigned int j=1; j<sequenceSize2+1; j++) { 
		for (unsigned int i=1; i<sequenceSize+1; i++) { 
			
			matrixE[i]=max(matrixH[i]-penalOpen,matrixE[i]-penalExtension);
			F=max(H_up-penalOpen, F-penalExtension);
			H_diag = matrixH[i-1] + m_matrixBlosum[databaseSequence[j-1]][querySequence[i-1]];
			H = max(H_diag, matrixE[i]);
			H = max(H, F);
			H = max(H, 0);
			
			if(H > scoreMax){
				scoreMax = H;
			}
			
			matrixH[i-1] = H_up;
			H_up= H;
		}
	
		matrixH[sequenceSize] = H;
		
		H_up = 0;
	}
	
	delete[] matrixH;
	delete[] matrixE;
	//delete[] databaseSequence;
	
	return scoreMax;
	
}

void Algorithm::calculate(Sequence * sequence, Database * database, int m_matrixBlosum[28][28]){
	
	// prendre les informations utiles de la database et de la query sequence
	int * querySequence = sequence->getSequence(); 
	int querySequenceSize = sequence->getSize(); 
	vector<uint32_t> sequenceOffsetTable = database->getSequenceOffsetTable(); 
	vector<uint32_t> headerOffsetTable = database->getHeaderOffsetTable();
	
	// ouvrir le sequence file
	ifstream sequenceFile;
	sequenceFile.open(database->getSequenceFileName());
	
	// initialiser les tableaux pour enregistrer les 50 meilleurs scores
	int results[50] = {0};
	int scores[50] = {0};
	
	// parcourir toutes les séquences de la database, calculer leur score et garder les 50 meilleurs scores
	if (sequenceFile.is_open()){
		for (unsigned int p = 0; p<sequenceOffsetTable.size()-1; p++){
			sequenceFile.seekg(sequenceOffsetTable[p]);
			int databaseSequenceSize = sequenceOffsetTable[p+1]-sequenceOffsetTable[p];
			int* databaseSequence = new int[databaseSequenceSize];
			for (unsigned int i=0; i<databaseSequenceSize; i++){
				char x;
				if (sequenceFile.get(x)){
					databaseSequence[i] = (x);
				}
			}
			
			int score = calculateScore(querySequence, querySequenceSize, databaseSequenceSize, m_matrixBlosum, databaseSequence);
			
			delete[] databaseSequence;
			
			// regarder si la séquence considérée a un meilleur alignement (meilleur score) que celles déjà retenues
			
			if (score > scores[0]){
				scores[0] = score;
				results[0] = p;
				int i = 0;
				while (scores[i] > scores[i+1]){
					int temp = scores[i];
					scores[i] = scores[i+1];
					scores[i+1] = temp;
					int temp2 = results[i];
					results[i] = results[i+1];
					results[i+1] = temp2;
					i++;
					if (i == 49){
						break;
					}
				}
			}
		}
	} else {
		cout << "You do not have the right sequence file name. The name of the file should be " << database->getSequenceFileName() << endl;
	}
	sequenceFile.close();
	// imprimer les résultats
	database->printSequenceName(results, scores);
}
