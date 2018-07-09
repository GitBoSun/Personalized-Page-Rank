#include <cstdio>
#include <cstdlib>
#include <vector>
#include <set>
#include <string.h>
#include <cmath>
#include <mpi.h>
#include <omp.h>


using std::vector;
using std::set;

#define MAX_ITER 10000
#define N 10000
#define MAX_OUT_NODE 20
#define STEP 100

typedef unsigned long ulong;

ulong *graph = NULL, begin_index[N + 1] = {}, count[N] = {}, total_count[N] = {}, old_total_count[N] = {};

const double alpha = 0.15;

double PR[N] = {};


void init_graph()
{
	srand(0);
	ulong size;
	begin_index[0] = 0;
	vector<ulong> g;
	set<ulong> out;
	for (ulong i = 0; i < N; ++i)
	{
		size = rand() % MAX_OUT_NODE;
		begin_index[i + 1] = begin_index[i] + size;
		while (out.size() < size)
		{
			out.emplace(rand() % N);
		}
		g.insert(g.end(), out.begin(), out.end());
	}
	graph = new ulong[g.size()];
	memcpy(graph, g.data(), g.size() * sizeof(ulong));
}


template <class T>
double mean(const T X[])
{
	double sum = 0.0;
	#pragma omp parallel for reduction(+:sum)
	for (int i = 0; i < N; ++i)
		sum += X[i];
	return sum / N;
}


template <class T>
double covariance(const T X[], const T Y[], const double x_mean, const double y_mean)
{
	double sum = 0.0;
	#pragma omp parallel for reduction(+:sum)
	for (int i = 0; i < N; i++)
		sum += (X[i] - x_mean) * (Y[i] - y_mean);
	return sum / (N - 1);
}


template <class T>
double variance(const T X[], const double x_mean)
{
	return covariance(X, X, x_mean, x_mean);
}


//返回样例的标准差
template <class T>
double standardDeviation(const T X[], const double x_mean)
{
	return sqrt(variance(X, x_mean));
}


//计算相关系数
template <class T>
double sampleCorrelationCoefficient(T X[], T Y[])
{
	double x_mean, y_mean, x_variance, y_variance, xy_covariance;
	#pragma omp parallel
	{
		#pragma omp task
		x_mean = mean(X);
		#pragma omp task
		y_mean = mean(Y);
		#pragma omp taskwait
		
		#pragma omp task
		x_variance = variance(X, x_mean);
		#pragma omp task
		y_variance = variance(Y, y_mean);
		#pragma omp task
		xy_covariance = covariance(X, Y, x_mean, y_mean);
		#pragma omp taskwait
	}
	return xy_covariance / sqrt(x_variance * y_variance);
}


int main()
{
	int MPI_Size, MPI_Rank;
	ulong current_vertex, out_size;
	clock_t start, end;
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &MPI_Size);
	MPI_Comm_rank(MPI_COMM_WORLD, &MPI_Rank);
	if (MPI_Rank == 0)
	{
		init_graph();
		printf("graph inited\n");
		printf("number of nodes: %lu\n", N);
		printf("number of edges: %lu\n", begin_index[N]);
		printf("average number of edges: %.2lf\n", double(begin_index[N]) / N);
	}
	start = clock();
	ulong p_t_iter = 0;
	MPI_Bcast(begin_index, N + 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
	if (MPI_Rank != 0)
	{
		graph = new ulong[begin_index[N]];
	}
	MPI_Bcast(graph, begin_index[N], MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
	int p_iter = MAX_ITER / MPI_Size;
	int p_step = STEP / MPI_Size;
	ulong t_iter;
	for (ulong source = 0; source < N; ++source)
	{
		if (begin_index[source] == begin_index[source + 1])
			continue;
		memset(count, 0, N * sizeof(ulong));
		memset(old_total_count, 0, N * sizeof(ulong));
		current_vertex = source;
		int iter;
		for (iter = 0; iter < p_iter; ++iter)
		{
			out_size = begin_index[current_vertex + 1] - begin_index[current_vertex];
			if (out_size < 1 || rand() < alpha * RAND_MAX)
				current_vertex = source; //因为是personalize的，所以有一定概率回到source，如果是pagerank，就从N个里面随机找一个
			else
			{
				current_vertex = graph[begin_index[current_vertex] + rand() % out_size];
			}
			++count[current_vertex];
			if ((iter + 1) % STEP == 0)
			{
				MPI_Reduce(count, total_count, N, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
				int flag = 0;
				if (MPI_Rank == 0)
				{
					if (iter + 1 > STEP)
					{
						double r = sampleCorrelationCoefficient(old_total_count, total_count);
						if (r > 0.999)
						{
							flag = 1;
						}
					}
					memcpy(old_total_count, total_count, N * sizeof(ulong));
				}
				MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
				if (flag != 0)
				{
					++iter;
					break;
				}
			}
		}
		p_t_iter += iter;
		MPI_Reduce(&iter, &t_iter, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(count, total_count, N, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		if (MPI_Rank == 0)
		{
			for (ulong i = 0; i < N; ++i)
				PR[i] += double(total_count[i]) / (t_iter * N);
		}
	}
	MPI_Reduce(&p_t_iter, &t_iter, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	delete[] graph;
	if (MPI_Rank == 0)
	{
		end = clock();
		printf("total iteration time: %lu\n", t_iter);
		printf("average iteration time: %.2lf\n", double(t_iter) / N);
		printf("total execution time: %.2lf ms\n", (1000.0 * (end - start)) / CLOCKS_PER_SEC);
	}
	MPI_Finalize();
}
