#include <mpi.h> 
#include <iostream>
#include <ctime>
#include <assert.h>
#include <cstdlib>
#include <cmath>
using namespace std;
void M_MIN(void *sb,void *rb,int n,MPI_Datatype t,MPI_Op op,int rt,MPI_Comm com) 
{
MPI_Status st;
int ProcNum;
MPI_Comm_size(com, &ProcNum);
for (int j = 0; j < ProcNum; j++) {
if (j != rt)
MPI_Recv(sb, n, t, j, 0, com, &st);
for (int i = 0; i < n; i++)
{
if (j == rt) {
((int *)rb)[i] = ((int *)sb)[i];
cout << "result in root " << ((int *)rb)[i] << endl;
}
else {
if (((int *)rb)[i] < ((int *)sb)[i])
((int *)rb)[i] = ((int *)sb)[i];
//int tmp = ((int *)recvbuf)[i];
cout << "result " << j << " " << ((int *)rb)[i] << "  last " << ((int *)sb)[i] << endl;
}
}
}
}
void M_MAX(void *sf, void *rf, int n, MPI_Datatype t, MPI_Op op, int rt, MPI_Comm com)
{
MPI_Status st;
int ProcNum;
MPI_Comm_size(com, &ProcNum);
for (int j = 0; j < ProcNum; j++) {
if (j != rt)
MPI_Recv(sf, n, t, j, 0, com, &st);
for (int i = 0; i < n; i++)
{
if (j == rt) {
((int *)rf)[i] = ((int *)sf)[i];
cout << "result in root " << ((int *)rf)[i] << endl;
}
else {
if (((int *)rf)[i] < ((int *)sf)[i])
((int *)rf)[i] = ((int *)sf)[i];
cout << "result " << j << " " << ((int *)rf)[i] << "  last " << ((int *)sf)[i] << endl;
}
}
}
}
void LXOR(void *sf,void *rf,int n,MPI_Datatype t,MPI_Op op,int rt,MPI_Comm com)
{
MPI_Status st;
int ProcNum;
MPI_Comm_size(com, &ProcNum);
for (int j = 0; j < ProcNum; j++) {
if (j != rt)
MPI_Recv(sf, n, t, j, 0, com, &st);
for (int i = 0; i < rt; i++)
{
if (j == rt) {
((int *)rf)[i] = ((int *)sf)[i];
cout << "result in root " << ((int *)rf)[i] << endl;
}
else {
((int *)rf)[i] = ((int *)rf)[i] != ((int *)sf)[i];
//int tmp = ((int *)recvbuf)[i];
cout << "result " << j << " " << ((int *)rf)[i] << "  last " << ((int *)sf)[i] << endl;
}
}
}
}
void BXOR(void *sf,void *rf,int n,MPI_Datatype t,MPI_Op op,int rt,MPI_Comm com)
{
MPI_Status st;
int ProcNum;
MPI_Comm_size(com, &ProcNum);
for (int j = 0; j < ProcNum; j++) {
if (j != rt)
MPI_Recv(sf, n, t, j, 0, com, &st);
for (int i = 0; i < n; i++)
{
if (j == rt) {
((int *)rf)[i] = ((int *)sf)[i];
cout << "result in root " << ((int *)rf)[i] << endl;
}
else {
((int *)rf)[i] = ((int *)rf)[i] ^ ((int *)sf)[i];
//int tmp = ((int *)recvbuf)[i];
cout << "result " << j << " " << ((int *)rf)[i] << "  last " << ((int *)sf)[i] << endl;
}
}
}
}
void BOR(void *sf,void *rf,int n,MPI_Datatype t,MPI_Op op,int rt,MPI_Comm com)
{
MPI_Status st;
int ProcNum;
MPI_Comm_size(com, &ProcNum);
for (int j = 0; j < ProcNum; j++) {
if (j != rt)
MPI_Recv(sf, n, t, j, 0, com, &st);
for (int i = 0; i < n; i++)
{
if (j == rt) {
((int *)rf)[i] = ((int *)sf)[i];
cout << "result in root " << ((int *)rf)[i] << endl;
}
else {
((int *)rf)[i] = ((int *)rf)[i] | ((int *)sf)[i];
//int tmp = ((int *)recvbuf)[i];
cout << "result " << j << " " << ((int *)rf)[i] << "  last " << ((int *)sf)[i] << endl;
}
}
}
}
void LOR(void *sf,void *rf,int n,MPI_Datatype t,MPI_Op op,int rt,MPI_Comm com)
{
MPI_Status st;
int ProcNum;
MPI_Comm_size(com, &ProcNum);
for (int j = 0; j < ProcNum; j++) {
if (j != rt)
MPI_Recv(sf, n, t, j, 0, com, &st);
for (int i = 0; i < n; i++)
{
if (j == rt) {
((int *)rf)[i] = ((int *)sf)[i];
cout << "result in root " << ((int *)rf)[i] << endl;
}
else {
((int *)rf)[i] = ((int *)rf)[i] || ((int *)sf)[i];
//int tmp = ((int *)recvbuf)[i];
cout << "result " << j << " " << ((int *)rf)[i] << "  last " << ((int *)sf)[i] <<endl;
}
}
}
}
void BAND(void *sf, void *rf, int n, MPI_Datatype t, MPI_Op op, int rt, MPI_Comm comm)
 {
MPI_Status st;
int ProcNum;
MPI_Comm_size(comm, &ProcNum);
for (int j = 0; j < ProcNum; j++) {
if (j != rt)
MPI_Recv(sf,n, t, j, 0, comm, &st);
for (int i = 0; i < n; i++)
{
if (j == rt) {
((int *)rf)[i] = ((int *)sf)[i];
cout << "result in root " << ((int *)rf)[i] << endl;
}
else {
((int *)rf)[i] = ((int *)rf)[i] & ((int *)sf)[i];

cout << "result " << j << " " << ((int *)rf)[i] << "  last " << ((int *)sf)[i] << endl;
}
}
}
}
void LAND(void *sf, void *rf, int n, MPI_Datatype t, MPI_Op op, int rt, MPI_Comm comm)
 {
MPI_Status st;
int ProcNum;
MPI_Comm_size(comm, &ProcNum);
for (int j = 0; j < ProcNum; j++) {
if (j != rt)
MPI_Recv(sf, n, t, j, 0, comm, &st);
for (int i = 0; i < n; i++)
{
if (j == rt) {
((int *)rf)[i] = ((int *)sf)[i];
cout << "result in root " << ((int *)rf)[i] << endl;
}
else {
((int *)rf)[i] = ((int *)rf)[i] && ((int *)sf)[i];

cout << "result " << j << " " << ((int *)rf)[i] << "  last " << ((int *)sf)[i] << endl;
}
}
}
}
void PROD(void *sf, void *rf, int n, MPI_Datatype t,MPI_Op op,int rt,MPI_Comm comm)
 {
MPI_Status st;
int ProcNum;
MPI_Comm_size(comm, &ProcNum);
for (int j = 0; j < ProcNum; j++) {
if (j != rt)
MPI_Recv(sf, n, t, j, 0, comm, &st);
for (int i = 0; i < n; i++)
{
if (j == rt) {
((int *)rf)[i] = ((int *)sf)[i];
cout << "result in root " << ((int *)rf)[i] << endl;
}
else {
((int *)rf)[i] = ((int *)rf)[i] * ((int *)sf)[i];
//int tmp = ((int *)recvbuf)[i];
cout << "result " << j << " " << ((int *)rf)[i] << "  last " << ((int *)sf)[i] << endl;
}
}
}
}
void MIN(void *sf, void *rf, int n, MPI_Datatype t,MPI_Op op,int rt,MPI_Comm comm)
 {
MPI_Status st;
int ProcNum;
MPI_Comm_size(comm, &ProcNum);
for (int j = 0; j < ProcNum; j++) {
if (j != rt)
MPI_Recv(sf, n, t, j, 0, comm, &st);
for (int i = 0; i < n; i++)
{
if (j == rt) {
((int *)rf)[i] = ((int *)sf)[i];
cout << "result in root " << ((int *)rf)[i] << endl;
}
else {
if (((int *)rf)[i] > ((int *)sf)[i])
((int *)rf)[i] = ((int *)sf)[i];
//int tmp = ((int *)rf)[i];
cout << "result " << j << " " << ((int *)rf)[i] << "  last " << ((int *)sf)[i] << endl;
}
}
}
}
void MAX(void *sf, void *rf, int n, MPI_Datatype t, MPI_Op op, int rt, MPI_Comm comm)
 {
MPI_Status st;
int ProcNum;
MPI_Comm_size(comm, &ProcNum);
for (int j = 0; j < ProcNum; j++) {
if (j != rt)
MPI_Recv(sf, n, t, j, 0, comm, &st);
for (int i = 0; i < n; i++)
{
if (j == rt) {
((int *)rf)[i] = ((int *)sf)[i];
cout << "result in root " << ((int *)rf)[i] << endl;
}
else {
if (((int *)rf)[i] < ((int *)sf)[i])
((int *)rf)[i] = ((int *)sf)[i];

cout << "result " << j << " " << ((int *)rf)[i] << "  last " << ((int *)sf)[i] << endl;
}
}
	}
}
void SUMM(void *sf, void *rf, int n, MPI_Datatype t, MPI_Op op,int rt,MPI_Comm comm)
 {
//std::cout<<"I am in 2 if" << std::endl;
MPI_Status st;
int ProcNum;
MPI_Comm_size(comm, &ProcNum);
for (int j = 0; j < ProcNum; j++) {
if (j != rt)
MPI_Recv(sf, n, t, j, 0, comm, &st);
for (int i = 0; i < n; i++)
{
if (j == rt) {
((int *)rf)[i] = ((int *)sf)[i];
//std::cout << "result in root " << ((int *)rf)[i] << endl;
}
else {
((int *)rf)[i] = ((int *)rf)[i] + ((int *)sf)[i];
//int tmp = ((int *)rf)[i];
}
}
}
}
void MY_MPI_SUMM_Tree(void *sendbuf, void *recvbuf, int count) {
//std::cout<<"I am in 2 if" << std::endl;
//MPI_Status st;
int ProcNum;
MPI_Comm_rank(MPI_COMM_WORLD, &ProcNum);
for (int i = 0; i < count; i++)
{
//if (j == root) {
//((int *)recvbuf)[i] = ((int *)sendbuf)[i];
//std::cout << "result in root " << ((int *)recvbuf)[i] << std::endl;
//}
//else {
((int *)recvbuf)[i] = ((int *)recvbuf)[i] + ((int *)sendbuf)[i];
//int tmp = ((int *)recvbuf)[i];
//}
}
}
int MY_MPI_Reduce(void *sf,void *rf,int n,MPI_Datatype t,MPI_Op op,int rt,MPI_Comm comm)
{
int ProcNum, ProcRank;
MPI_Comm_size(comm, &ProcNum);
MPI_Comm_rank(comm, &ProcRank);
if (ProcRank != rt) {
//std::cout << "I am in not root" << std::endl;
MPI_Send(sf,n, t, rt, 0, comm);
}
else {
if (ProcRank == rt)
{
//std::cout << "I am in root " << ProcNum << std::endl;
//	MPI_Status st;
//int j = 0;
if (t == MPI_INT) {
if (op == MPI_SUM) {
SUMM(sf, rf, n, t, op, rt, comm);
return 0;
}
if (op == MPI_MAX) {
M_MAX(sf, rf, n, t, op, rt, comm);
return 0;
}
if (op == MPI_MIN) {
M_MIN(sf, rf, n, t, op, rt, comm);
return 0;
}
if (op == MPI_PROD) {
PROD(sf, rf, n, t, op, rt, comm);
return 0;
}
if (op == MPI_LAND) {
LAND(sf, rf, n, t, op, rt, comm);
return 0;
}
if (op == MPI_LOR) {
LOR(sf, rf, n, t, op, rt, comm);
return 0;
}
if (op == MPI_BAND) {
BAND(sf, rf, n, t, op, rt, comm);
return 0;
}
if (op == MPI_BOR) {
BOR(sf, rf, n, t, op, rt, comm);
return 0;
}
if (op == MPI_BXOR) {
BXOR(sf, rf, n, t, op, rt, comm);
return 0;
}
if (op == MPI_LXOR) {
LXOR(sf, rf, n, t, op, rt, comm);
}
if (op == MPI_MAXLOC)
return -1;
//MY_MPI_MAXLOC(sendbuf, recvbuf, count, type, op, root, comm);
if (op == MPI_MINLOC)
return -1;
}
}
}
 return -10;
}

void rec(void *sf,void *rf,int n,MPI_Datatype t,MPI_Op op,int rt,MPI_Comm comm,int* mas,int size,int h)
 {
int curentSize, procNum, rank;
int* curentMass;
int* countMassElement = new int[n];
MPI_Comm_size(comm, &procNum);
MPI_Comm_rank(comm, &rank);
MPI_Status st;
if (size % 2) {
curentSize = size / 2 + 1;
curentMass = new int[curentSize];
}
else {
curentSize = size / 2;
curentMass = new int[curentSize];
}
for (int i = 0; i < curentSize; i++) {
curentMass[i] = -1;
}
h++;
//int t = 0;
for (int i = 0, t = 0; i < size; i++) {
if (mas[i] == rank) {
if (!(size % 2)) {
if (!(i % 2)) {
MPI_Send(sf, n, t, mas[i + 1], 0, comm);
for (int q = 0; q < n; q++)
{
countMassElement[q] = ((int *)(rf))[q];
}
}
else {
MPI_Recv(rf, n, t, mas[i - 1], 0, comm, &st);
MY_MPI_SUMM_Tree(sf, rf, n);
sf = rf;
}
}
else {
if ((rank != procNum - 1)) {
if (!(i % 2)) {
MPI_Send(sf, n, t, mas[i + 1], 0, comm);
for (int q = 0; q < n; q++)
{
countMassElement[q] = ((int *)(rf))[q];
}
}
else {
for (int r = 0; r < n; r++) {
((int *)rf)[i] = ((int *)sf)[i];
}
MPI_Recv(sf, n, t, mas[i - 1], 0, comm, &st);
MY_MPI_SUMM_Tree(sf, rf, n);
for (int r = 0; r < n; r++) {
((int *)sf)[i] = ((int *)rf)[i];
}
}
}
else {
curentMass[t] = i;
t++;
}
}
}
if (!(size % 2)) {
if (i % 2) {
curentMass[t] = i;
t++;
}
}
}
if (curentSize != 1) {
rec(sf, countMassElement, n, t, op, rt, comm, curentMass, curentSize, h);
}
else {
if (curentMass[0] != rt) {
if (rank == curentMass[0])
{
MPI_Send(sf, n, t, rt, 0, comm);
}
if (rank == rt) {
MPI_Recv(rf, n, t, curentMass[0], 0, comm, &st);
}
}
else {
rf = sf;
}
}
}
int Tree(void *sf, void *rf, int n, MPI_Datatype t, MPI_Op op, int rt, MPI_Comm com)
{
int ProcNum, ProcRank;
MPI_Comm_size(com, &ProcNum);
MPI_Comm_rank(com, &ProcRank);
int * massProcRankSend = new int[ProcNum];
for (int i = 0; i < ProcNum; i++) {
massProcRankSend[i] = i;
}
	
if (ProcRank == rt) {
}
return 0;
}
int main(int argc, char* argv[]) {
int n = atoi(argv[1]);
int *mas = new int[n];
int *mas_r = new int[n];
int ProcRank;
double Time_my1 = 0, Time_my2 = 0,Time_tree1=0, Time_tree2=0, Time_MPI1 = 0, Time_MPI2 = 0;
MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
srand(static_cast<int>(time(NULL)));
for (int i = 0; i < n; i++) {
mas[i] = 10 + std::rand() % 1000 + ProcRank;
}
if (ProcRank == 0)
{
Time_my1 = MPI_Wtime();
}
MY_MPI_Reduce(mas, mas_r, n, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
if (ProcRank == 0)
{
Time_my2 = MPI_Wtime();
Time_tree1 = MPI_Wtime();
}
Tree(mas, mas_r, n, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
if (ProcRank == 0)
{
Time_tree2 = MPI_Wtime();
//std::cout << " mass resualt ";
for (int i = 0; i < n; i++) {
//std::cout << mas_r[i] << " ";
}
//cout << std::endl;
Time_MPI1 = MPI_Wtime();
}
	
MPI_Reduce(mas, mas_r, n, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
if (ProcRank == 0)
{
	
Time_MPI2 = MPI_Wtime();
cout << "reduce realisation MPI: t =  " << Time_MPI2 - Time_MPI1 << endl;
cout << "reduce realisation my tree: t =  " << Time_tree2 - Time_tree1 << endl;
cout << "reduce realisation my: t =  " << Time_my2 - Time_my1 << endl;
}
MPI_Finalize();
return 0;
}