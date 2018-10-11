#include <iostream>
#include <mpi.h>
#include <string>
#include <ctime>
using namespace std;

int main(int argc, char *argv[])
{
	string s = "hngb hjnhgb jhbjgb hjbbj hhh hhhhjbvhg  hjbjbjgb j j jjjjj";
	int len = s.length();

	int ProcNum, ProcRank;
	double stime = 0.0;
	double etime = 0.0;
	int nword = 0;
	int nresword = 0;
	int residue = 0;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	if (ProcRank == 0)
	{
		stime = MPI_Wtime();
		for (int i = 0; i < len; i++)
			if (s[i] == ' ')
				nword++;
		etime = MPI_Wtime();
		cout << "1 process" << endl;
		cout << "Number words: " << nword << endl;
		cout << "Time: " << etime - stime << " sec" << endl;
		nword = 0;
		stime = MPI_Wtime();
	}

	MPI_Bcast(s, len, MPI_CHAR, 0, MPI_COMM_WORLD);

	for (int i = len / ProcNum * ProcRank; i < len / ProcNum + len / ProcNum * ProcRank; i++)
		if (s[i] = ' ')
			nword++;

	MPI_Reduce(&nword, &nresword, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if (ProcRank == 0)
	{
		residue = len % ProcNum;
		if (residue)
			for (int i = len - residue; i < len; i++)
				if (s[i] == ' ')
					nresword++;
		etime = MPI_Wtime();
		cout << ProcNum << " process" << endl;
		cout << "Number words: " << nresword << endl;
		cout << "Time: " << etime - stime << " sec" << endl;
	}

	MPI_Finalize();
	return 0;
}