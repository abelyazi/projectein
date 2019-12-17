#include "Algorithm.h"
#define max(x,y) ( x < y ? y : x ) //

Algorithm::Algorithm(int R, int Q){
	penalExtension = R; // R
	penalOpen = Q; // Q
}

int Algorithm::calculateScore(vector<int> querySequence, int sequenceSize, int sequenceSize2, int m_matrixBlosum[28][28], int * databaseSequence){
	
	int* matrixH = new int[sequenceSize+1];
	int* matrixE = new int[sequenceSize+1];
	int F = 0;
	int Hhaut = 0;
	int Hdiag = 0;
	int Hmaintenant = 0;
	int scoreMax =0;
	
	for (unsigned int i=0; i< sequenceSize+1;i++)
		{
		matrixE[i] = 0;
		matrixH[i] = 0;
		
		}

	for (unsigned int j=1; j<sequenceSize2+1; j++) { 
		for (unsigned int i=1; i<sequenceSize+1; i++) { 
			matrixE[i]=max(matrixH[i]-penalOpen,matrixE[i]-penalExtension);
			F=max(Hhaut-penalOpen, F-penalExtension);
			Hdiag = matrixH[i-1] + m_matrixBlosum[databaseSequence[j-1]][querySequence[i-1]];
			Hmaintenant = max(Hdiag, matrixE[i]);
			Hmaintenant = max(Hmaintenant, F);
			Hmaintenant = max(Hmaintenant, 0);

			if(Hmaintenant > scoreMax){
				scoreMax = Hmaintenant;
			}
			matrixH[i-1] = Hhaut;
			Hhaut= Hmaintenant;
		}
	
		matrixH[sequenceSize] = Hmaintenant;
		
		Hhaut = 0;
	}
	
	delete[] matrixH;
	delete[] matrixE;
	delete[] databaseSequence;
	
	
	
	
	
	
	
	
	/*
	int* matrixH = new int[sequenceSize+1];
	int* matrixE = new int[sequenceSize+1];
	int F = 0;
	int Hhaut = 0;
	int Hdiag = 0;
	int Hmaintenant = 0;
	int scoreMax = 0;

	for (unsigned int i=0; i< sequenceSize+1;i++) {
		matrixE[i] = 0;
		matrixH[i] = 0;
	}

	for (unsigned int j=1; j<sequenceSize2+1; j++) {
		for (unsigned int i=1; i<sequenceSize+1; i++) {  
			matrixE[i]=max(matrixH[i]-penalOpen,matrixE[i]-penalExtension);
			F=max(Hhaut-penalOpen, F-penalExtension);
			Hdiag = matrixH[i-1] + m_matrixBlosum[databaseSequence[j-1]][querySequence[i-1]];
			Hmaintenant = max(Hdiag, matrixE[i]);
			Hmaintenant = max(Hmaintenant, F);
			Hmaintenant = max(Hmaintenant, 0);

			if(Hmaintenant > scoreMax){
				scoreMax = Hmaintenant;
			}
			matrixH[i-1] = Hhaut;
			Hhaut= Hmaintenant;
		}
	
		matrixH[sequenceSize] = Hmaintenant;
		
		Hhaut = 0;
	}
	delete[] matrixH;
	delete[] matrixE;
	delete[] databaseSequence;//
	
	/*
	// impression de la séquence trouvée dans le sequence file
	cout << "Séquence trouvée dans la database = ";
	for (int n:foundSequence){
		cout << n;
	} cout << endl;*/
	
	/*
	// créer les matrice H, E et F
	int* matrixH = new int[sequenceSize*sequenceSize2];
	int* matrixE = new int[sequenceSize*sequenceSize2];
	int* matrixF = new int[sequenceSize*sequenceSize2];
	for (unsigned int i=0; i< sequenceSize*sequenceSize2;i++) {
		matrixE[i] = 0;
		matrixF[i] = 0;
		matrixH[i] = 0;
		}		
	// calculer les scores dans la matrice H avec les matrices E et F
	for (unsigned int i=1; i<sequenceSize; i++) {
		for (unsigned int j=1; j<sequenceSize2; j++) {
			//Eij
			if(matrixE[i*sequenceSize2 + j-1] - penalExtension < matrixH[i*sequenceSize2 + j-1] - penalOpen ) {		
				matrixE[i*sequenceSize2 + j] = matrixH[i*sequenceSize2 + j-1] - penalOpen;
			}
			else {
				matrixE[i*sequenceSize2 + j] = matrixE[i*sequenceSize2 + j-1] - penalExtension;
			}
			//FIJ 
			if(matrixF[(i-1)*sequenceSize2 + j] - penalExtension < matrixH[(i-1)*sequenceSize2 + j] - penalOpen ) {		
				matrixF[i*sequenceSize2 + j] = matrixH[(i-1)*sequenceSize2 + j] - penalOpen;
			}
			else {
				matrixF[i*sequenceSize2 + j] = matrixF[(i-1)*sequenceSize2 + j] - penalExtension;
			}
			//Hij
			matrixH[i*sequenceSize2 + j] = matrixH[(i-1)*sequenceSize2+j-1] + m_matrixBlosum[querySequence[i-1]][databaseSequence[j-1]];
			if(matrixH[i*sequenceSize2 + j] < (matrixE[i*sequenceSize2 + j])) {
				matrixH[i*sequenceSize2 + j] = (matrixE[i*sequenceSize2 + j]);
			}
			if(matrixH[i*sequenceSize2 + j] < (matrixF[i*sequenceSize2 + j])) {
				matrixH[i*sequenceSize2 + j] = (matrixF[i*sequenceSize2 + j]);
			}
			if(matrixH[i*sequenceSize2 + j] < 0) {
				matrixH[i*sequenceSize2 + j] = 0;
			}
			
			// chercher le score maximum (= meilleur alignement possible) dans la matrice H
			if (matrixH[i*sequenceSize2 + j] > scoreMax){
				scoreMax = matrixH[i*sequenceSize2 + j];
			};
		}	
	}
	
	delete matrixH;
	delete matrixE;
	delete matrixF;
	*/
	
	return scoreMax;
	
}

void Algorithm::calculate(Sequence * sequence, Database * database, int m_matrixBlosum[28][28]){
	
	// ouvrir le sequence file
	ifstream sequenceFile;
	sequenceFile.open("./newE.fasta.psq");
	
	// prendre les informations utiles de la database et de la query sequence
	vector<int> querySequence = sequence->getSequence(); 
	int querySequenceSize = sequence->getSize(); 
	vector<uint32_t> sequenceOffsetTable = database->getSequenceOffsetTable(); 
	vector<uint32_t> headerOffsetTable = database->getHeaderOffsetTable(); 
	
	int results[5] = {0};
	int scores[5] = {0};
	
	// parcourir toutes les séquences de la database, calculer leur score et garder les 5 meilleurs scores
	for (unsigned int p = 0; p<sequenceOffsetTable.size()-1; p++){
		sequenceFile.seekg(sequenceOffsetTable[p]);
		int databaseSequenceSize = sequenceOffsetTable[p+1]-sequenceOffsetTable[p];
		//int databaseSequence[databaseSequenceSize] = {0};//
		int* databaseSequence = new int[databaseSequenceSize]; //
		for (unsigned int i=0; i<databaseSequenceSize; i++){
			char x;
			if (sequenceFile.get(x)){
				databaseSequence[i] = (x);
			}
		} 
		int score = calculateScore(querySequence, querySequenceSize, databaseSequenceSize, m_matrixBlosum, databaseSequence);
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
				if (i == 4){
					break;
				}
			}
		}
	}
	sequenceFile.close();
	database->printSequenceName(results, scores);
}
