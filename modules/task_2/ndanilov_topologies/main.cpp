#include <assert.h>
#include <time.h>
#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <string>

#include "host.h"

int main(int argc, char* argv[]) {
  /* Mpi variables */
  int status,
      rank,
      size,
      proc_num;
  int source, destination;
  std::string msg;

  if (argc > 3) {
    source = atoi(argv[1]);
    destination =  atoi(argv[2]);
    msg = std::string(argv[3]);
  } else {
    source = -1;
    destination = -1;
    msg = "Hi";
  }

  assert(source >= 0);
  assert(destination >= 0);

  /* Mpi init block */
  status = MPI_Init(&argc, &argv);

  status = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  status = MPI_Comm_size(MPI_COMM_WORLD, &size);

  ring_host *test_host = new ring_host(rank, size);

  if (source == rank) {
    assert(msg != std::string() || msg.size() < 255);
    test_host->generate_packet(source, destination, msg);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  test_host->xmit();

  return 0;
}
