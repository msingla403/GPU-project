#include<bits/stdc++.h>
#include<cuda.h>
#include <thrust/sort.h>
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

#define PANEL_SIZE 3
#define DENSE_THRESHOLD 2

// __global__ void SPMM(int * tile_row_ptr, int * panel_ptr, int * col_val, int * col_idx){

// 	int row_panel_id = blockIdx.x;
// 	int row_id = threadIdx.x/32;
// 	int thread_no = threadIdx.x%32;

// 	int num_tiles = panel_ptr[row_panel_id+1] - panel_ptr[row_panel_id];

// 	int ptr = panel_ptr[row_panel_id]*PANEL_SIZE + row_id*num_tiles;

// 	for(int i=0;i<num_tiles;++i){

// 		int low = tile_row_ptr[ptr+i];
// 		int high = tile_row_ptr[ptr+i+1];

// 		if(high>low){
// 			int j=low;
// 			O[row_id][thread_no] += col_val[j] * D[col_idx[j]][thread_no];
// 		}
// 	}
// }

__global__ void find_dense(int *col_ptr, int* row_idx, int *isdense, int nr, int nc)
{	
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	if(idx < (nr/PANEL_SIZE)*nc)
	{
		int panel_id = idx/nc;
		int col_id = idx%nc;
		
		int counter = 0;
		for(int i=col_ptr[col_id]; i<col_ptr[col_id+1]; i++)
		{
			if(row_idx[i] >= panel_id*PANEL_SIZE && row_idx[i] < (panel_id+1)*PANEL_SIZE)
			{	
				// printf("%d %d %d\n", idx, col_id, row_idx[i]);
				counter++;
			}
		}
		// __syncthreads();
		// printf("%d %d %d\n", idx, counter, panel_id);
		if(counter >= DENSE_THRESHOLD)
			isdense[idx] = 1;
		else
			isdense[idx] = 0;
	}

}

int main(int argc, char** argv)
{
	char* inputfilename = argv[1];
	FILE *fp;
	fp = fopen(inputfilename, "r");
	
	int nr, nc, ne;
	fscanf(fp, "%d %d %d", &nr, &nc, &ne);

	
	int r, c;
	int rows[ne];
	int cols[ne];
	for(int i=0; i<ne; i++)
	{
		fscanf(fp, "%d %d", &r, &c);
		
		rows[i] = r;
		cols[i] = c; 
	}
	
	
	// // create column wise CSR
	thrust::sort_by_key(cols, cols+ne, rows);
	vi col_ptr(nc+1, 0);
	vi row_idx(ne, 0);
	for(int i=0; i<ne; i++)
	{
		// cout << cols[i] << " " << rows[i] << endl;
		col_ptr[cols[i]]++;
		row_idx[i] = rows[i]-1;
	}
	
	for(int i=0; i<nc; i++)
	col_ptr[i+1] += col_ptr[i];
	
	// for(int i=0; i<=nc; i++)
	// 	cout << col_ptr[i] << " ";
	// cout <<endl;
	
	// for(int i=0; i<ne; i++)
	// 	cout << row_idx[i] << " ";
	// cout << endl;
	
	// // Create row wise CSR
	thrust::sort_by_key(rows, rows+ne, cols);
	vi row_ptr(nr+1, 0);
	vi col_idx(ne, 0);
	
	// int col_val[ne];
	for(int i=0; i<ne; i++)
	{
		// cout << rows[i] << " " << cols[i] << endl;
		row_ptr[rows[i]]++;
		col_idx[i] = cols[i]-1;
	}
	
	for(int i=0; i<nr; i++)
		row_ptr[i+1] += row_ptr[i];


	// find dense tiles now
	int num_panels = nr/PANEL_SIZE;
	int thr = num_panels*nc;

	int *dcol_ptr;
	int *drow_idx;
	int *is_dense;
	cudaMalloc(&dcol_ptr, (nc+1)*sizeof(int));
	cudaMalloc(&drow_idx, ne*sizeof(int));
	cudaMalloc(&is_dense, thr*sizeof(int));

	cudaMemcpy(dcol_ptr, &col_ptr[0], (nc+1)*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(drow_idx, &row_idx[0], ne*sizeof(int), cudaMemcpyHostToDevice);

	find_dense <<<(thr + 1023)/1024, 1024>>>(dcol_ptr, drow_idx, is_dense, nr, nc);

	int isdense[thr];
	cudaMemcpy(isdense, is_dense, thr*sizeof(int), cudaMemcpyDeviceToHost);
	cudaFree(dcol_ptr);
	cudaFree(drow_idx);
	cudaFree(is_dense);

	// for(int i=0; i<num_panels; i++)
	// {
	// 	for(int j=0; j<nc; j++)
	// 		cout << isdense[i*nc + j] << " ";
	// 	cout << endl;
	// }
	

	// // Reorder CSR so that dense elements are before the sparse ones
	vi panel_ptr(num_panels+1, 0);
		vi tile_row_ptr;
	for(int panel_id=0; panel_id<num_panels; ++panel_id)
	{
		set <int> densecols;
		for(int j=0; j<nc; j++)
		{
			if(isdense[panel_id*nc + j])
				densecols.insert(j);
		}
		panel_ptr[panel_id+1] = densecols.size()+1; // one sparse panel

		for(int i = panel_id*PANEL_SIZE; i<(panel_id+1)*PANEL_SIZE; ++i)
		{
			if(i >= nr)
				break;

			vi temp1;
			vi temp2;
			for(int k=row_ptr[i]; k<row_ptr[i+1]; ++k)
			{
				if(densecols.find( col_idx[k] ) == densecols.end() )
					temp2.push_back(col_idx[k]);
				else
					temp1.push_back(col_idx[k]);
			}

			int counter = 0;
			for(int k=row_ptr[i]; k<row_ptr[i+1]; ++k)
			{
				if(counter < temp1.size())
				{
					col_idx[k] = temp1[counter];
				}
				else
				{
					col_idx[k] = temp2[counter-temp1.size()];
				}
				++counter;
			}
		}
		// densecols.clear();
	}

	for(int i=0; i<num_panels; i++)
		panel_ptr[i+1] += panel_ptr[i];

	// for(int i=0; i<=nr; i++)
	// 	cout << row_ptr[i] << " ";
	// cout <<endl;
	
	// for(int i=0; i<ne; i++)
	// 	cout << col_idx[i] << " ";
	// cout << endl;

	// for(int i=0; i<= num_panels; i++)
	// 	cout << panel_ptr[i] << " ";
	// cout << endl;

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



