#ifndef PACKET_H
#define PACKET_H

#define MAX_DATA_SIZE 255

struct packet {
	int src = -1;
	int dst = -1;
	char data[MAX_DATA_SIZE];
};

#endif
