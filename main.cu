#include<bits/stdc++.h>
#include<cuda.h>
using namespace std;

#define pairi pair<int,int>
#define ve vector
#define vi vector<int>
#define f first
#define s second
#define t third



__device__ __host__ int hashFn(int* data, int bsize)
{
	int res = bsize;
	for(int i=0; i<bsize; i++)
	{
		res ^= data[i] + 0x9e3779b9 + (res<<6) + (res>>2);
	}
	return abs(res);
}

__global__ void getSig(int *rowptr, int *colidx, int* perms, int* sigs, int siglen, int n)
{
	int idx =  blockIdx.x*blockDim.x + threadIdx.x;
	
	// if(idx == 0)
	// {
	// 	for(int i=0; i<n*siglen; i++)
	// 		printf("%d ", perms[i]);
	// 	printf("\n");
	// }

	if(idx <n)
	{
		for(int k=0; k<siglen; k++)
		{	
			int smallest = INT_MAX;
			for(int j=rowptr[idx]; j<rowptr[idx+1]; j++)
			{

				smallest = min(smallest, perms[k*n + colidx[j]]);
			}
			sigs[idx*siglen + k] = smallest;
		}
		// for(int i=0; i<siglen; i++)
		// {
		// 	printf("%d %d\n", idx, sigs[idx*siglen + i]);
		// }
	}  	
}

__global__ void getBuckets(int *sigs, int *res, int n, int siglen, int bsize, int numbuckets)
{
	int idx =  blockIdx.x*blockDim.x + threadIdx.x;

	if(idx < n)
	{
		int num_bands = siglen/bsize;
		for(int i=0; i<num_bands; i++)
		{
			int bkt = hashFn(&sigs[idx*siglen + i*bsize], bsize);
			res[idx*num_bands + i] = bkt%numbuckets;
		}
	}
}

set<pairi> LSH(vi &rowptr, vi &colidx, int siglen, int bsize, int numbuckets){
	int n = rowptr.size() - 1;

	int hperms[n*siglen];
	for(int k=0; k<siglen; k++)
	{
		vi perm(n);
		for(int i=0; i<n; i++)
		perm[i] = i;
		
		random_shuffle(perm.begin(), perm.end());
		copy(perm.begin(), perm.end(), &hperms[n*k]);		
	}

	int *drowptr;
	int *dcolidx;
	int *dperms;
	int *dsigs;
	cudaMalloc(&drowptr, rowptr.size()*sizeof(int));
	cudaMalloc(&dcolidx, colidx.size()*sizeof(int));
	cudaMalloc(&dperms, n*siglen*sizeof(int));
	cudaMalloc(&dsigs, n*siglen*sizeof(int));

	cudaMemcpy(drowptr, &rowptr[0], rowptr.size()*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dcolidx, &colidx[0], colidx.size()*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dperms, hperms, n*siglen*sizeof(int), cudaMemcpyHostToDevice);
	
	getSig<<< (n+1023/1024), 1024>>> (drowptr, dcolidx, dperms, dsigs, siglen, n);

	int sigs[n*siglen];
	cudaMemcpy(sigs, dsigs, n*siglen*sizeof(int), cudaMemcpyDeviceToHost);
	cudaFree(drowptr);
	cudaFree(dcolidx);
	cudaFree(dperms);
	// cudaFree(dsigs);

	// for(int i=0; i<n; i++)
	// {
	// 	for(int k=0; k<siglen; k++)
	// 		cout << sigs[i*siglen + k] << " ";
	// 	cout << endl;
	// }
	// cout << endl;
	
	int num_bands = siglen/bsize;
	int *dbucks;
	cudaMalloc(&dbucks, n*num_bands*sizeof(int));

	getBuckets<<<(n+1023)/1024, 1024>>>(dsigs, dbucks, n, siglen, bsize, numbuckets);
	int hbucks[n*num_bands];
	cudaMemcpy(hbucks, dbucks, n*num_bands*sizeof(int), cudaMemcpyDeviceToHost);
	
	vector<set<int>> buckets(numbuckets);
	for(int i=0; i<n; i++)
	{
		for(int j=0; j<num_bands; j++)
		{
			int idx = hbucks[i*num_bands + j];
			buckets[idx].insert(i);
		}
	}

	set<pairi> result;
	for(auto s: buckets)
	{
		vi temp(s.begin(), s.end());
		for(int i=0; i<temp.size(); i++)
		{
			for(int j=i+1; j<temp.size(); j++)
			{
				result.insert(make_pair(temp[i], temp[j]));
			}
		}
	}
	return result;
}

// struct posi{
// 	int r, c;
// };

int main(int argc, char** argv)
{
	char* inputfilename = argv[1];
	FILE *fp;
	fp = fopen(inputfilename, "r");
	
	int nr, nc, ne;
	fscanf(fp, "%d %d %d", &nr, &nc, &ne);

	int row_ptr[nr+1];
	for(int i=0; i<=nr; i++)
		row_ptr[i] = 0;

	int col_idx[ne];
	int col_val[ne];

	int r, c;
	
	for(int i=0; i<ne; i++)
	{
		fscanf(fp, "%d %d", &r, &c);
		
		row_ptr[r]++;
		col_idx[i] = c-1;
		col_val[i] = 1;		
	}
	for(int i=0; i<nr; i++)
		row_ptr[i+1] += row_ptr[i];
	
	for(int i=0; i<=nr; i++)
	{
		cout << row_ptr[i] << " ";
	}
	cout <<endl;

	for(int i=0; i<ne; i++)
		cout << col_idx[i] << " ";
	cout << endl;

	// int n = 6;
	// int m = 6;
	// vi rowptr{0,2,5,7,8,11,13};
	// vi colidx{0,4,1,3,5,2,4,1,0,3,4,2,5};
	
	// set<pairi> candidates = LSH(rowptr, colidx, 6, 2, 6);

	// for(auto i:candidates)
	// {
	// 	cout << i.f << " " << i.s << endl;
	// }
}



