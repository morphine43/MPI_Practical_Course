#ifndef HOST_H
#define HOST_H

#include <string>

#include "packet.h"

class host {
	public:

		int id = 0;
		packet pkt;
		MPI_Datatype type;
		int proc_number;

		host(int rank, int size);
		void generate_packet(int, int, std::string);
		virtual ~host();
		virtual int who_next() const = 0;
		virtual void xmit() = 0;
};

class ring_host : public host {
	public:
		int checked = 0;
		ring_host(int rank, int size);
		~ring_host();
		int who_next() const final;
		void xmit() final;
};
#endif