/* Copyright Nikita Danilov */
#ifndef MODULES_TASK_2_NDANILOV_TOPOLOGIES_INCLUDE_PACKET_H_
#define MODULES_TASK_2_NDANILOV_TOPOLOGIES_INCLUDE_PACKET_H_

#define MAX_DATA_SIZE 255

struct packet {
  int src = -1;
  int dst = -1;
  char data[MAX_DATA_SIZE];
};

#endif  // MODULES_TASK_2_NDANILOV_TOPOLOGIES_INCLUDE_PACKET_H_
