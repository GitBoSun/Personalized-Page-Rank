#include <cstdio>
#include <cstdlib>
#include <vector>
#include <set>
#include <string.h>
#include <cmath>
#include <ctime>
#include <iostream>


using std::vector;
using std::set;

#define MAX_ITER 10000
#define N 10000
#define MAX_OUT_NODE 20
#define STEP 100

typedef unsigned long ulong;

ulong *graph = NULL, begin_index[N + 1] = {}, count[N] = {}, old_count[N] = {};

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
	for (int i = 0; i < N; ++i)
		sum += X[i];
	return sum / N;
}


template <class T>
double covariance(const T X[], const T Y[], const double x_mean, const double y_mean)
{
	double sum = 0.0;
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
	double x_mean = mean(X);
	double y_mean = mean(Y);
	double x_variance = variance(X, x_mean);
	double y_variance = variance(Y, y_mean);
	double xy_covariance = covariance(X, Y, x_mean, y_mean);
	return xy_covariance / sqrt(x_variance * y_variance);
}


int main()
{
	ulong current_vertex, out_size;
	clock_t start, end;
	init_graph();
	printf("graph inited\n");
	printf("number of nodes: %lu\n", N);
	printf("number of edges: %lu\n", begin_index[N]);
	printf("average number of edges: %.2lf\n", double(begin_index[N]) / N);
	start = clock();
	ulong total_iteration_time = 0;
	for (ulong source = 0; source < N; ++source)
	{
		if (begin_index[source] == begin_index[source + 1])
			continue;
		memset(count, 0, N * sizeof(ulong));
		memset(old_count, 0, N * sizeof(ulong));
		current_vertex = source;
		int iter;
		for (iter = 0; iter < MAX_ITER; ++iter)
		{
			out_size = begin_index[current_vertex + 1] - begin_index[current_vertex];
			if (out_size == 0 || rand() < alpha * RAND_MAX)
				current_vertex = source; // 因为是 personalize 的，所以有一定概率回到 source，如果是 pagerank，就从N个里面随机找一个
			else
			{
				current_vertex = graph[begin_index[current_vertex] + rand() % out_size];
			}
			++count[current_vertex];
			if ((iter + 1) % STEP == 0)
			{
				if (iter + 1 > STEP)
				{
					double r = sampleCorrelationCoefficient(old_count, count);
					if (r > 0.999)
					{
						++iter;
						break;
					}
				}
				memcpy(old_count, count, N * sizeof(ulong));
			}
		}
		total_iteration_time += iter;
		for (ulong i = 0; i < N; ++i)
			PR[i] += double(count[i]) / (iter * N);
	}
	delete[] graph;
	end = clock();
	printf("total iteration time: %lu\n", total_iteration_time);
	printf("average iteration time: %.2lf\n", double(total_iteration_time) / N);
	printf("total execution time: %.2lf ms\n", (1000.0 * (end - start)) / CLOCKS_PER_SEC);
	return 0;
}
