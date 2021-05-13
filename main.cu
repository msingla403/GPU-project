#include<iostream>
#include<stdlib.h>
#include<set>
#include<map>
#include<vector>
#include<algorithm>
#include<cuda.h>
#include <thrust/sort.h>
using namespace std;

#define pairi pair<int,int>
#define ve vector
#define vi vector<int>
#define f first
#define s second
#define t third

#define PANEL_SIZE 21
#define DENSE_THRESHOLD 5

__device__ __host__ int hashFn(int *data, int bsize) {
	int res = bsize;
	for (int i = 0; i < bsize; i++) {
		res ^= data[i] + 0x9e3779b9 + (res << 6) + (res >> 2);
	}
	return abs(res);
}

__global__ void getSig(int *rowptr, int *colidx, int *perms, int *sigs, int siglen, int n) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	// if(idx == 0)
	// {
	// 	for(int i=0; i<n*siglen; i++)
	// 		printf("%d ", perms[i]);
	// 	printf("\n");
	// }

	if (idx < n) {
		for (int k = 0; k < siglen; k++) {
			int smallest = INT_MAX;
			for (int j = rowptr[idx]; j < rowptr[idx + 1]; j++) {

				smallest = min(smallest, perms[k * n + colidx[j]]);
			}
			sigs[idx * siglen + k] = smallest;
		}
		// for(int i=0; i<siglen; i++)
		// {
		// 	printf("%d %d\n", idx, sigs[idx*siglen + i]);
		// }
	}
}

__global__ void getBuckets(int *sigs, int *res, int n, int siglen, int bsize, int numbuckets) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < n) {
		int num_bands = siglen / bsize;
		for (int i = 0; i < num_bands; i++) {
			int bkt = hashFn(&sigs[idx * siglen + i * bsize], bsize);
			res[idx * num_bands + i] = bkt % numbuckets;
		}
	}
}

set<pairi > LSH(vi &rowptr, vi &colidx, int siglen, int bsize, int numbuckets) {
	int n = rowptr.size() - 1;

	int hperms[n * siglen];
	for (int k = 0; k < siglen; k++) {
		vi perm(n);
		for (int i = 0; i < n; i++)
			perm[i] = i;

		random_shuffle(perm.begin(), perm.end());
		copy(perm.begin(), perm.end(), &hperms[n * k]);
	}

	int *drowptr;
	int *dcolidx;
	int *dperms;
	int *dsigs;
	cudaMalloc(&drowptr, rowptr.size() * sizeof(int));
	cudaMalloc(&dcolidx, colidx.size() * sizeof(int));
	cudaMalloc(&dperms, n * siglen * sizeof(int));
	cudaMalloc(&dsigs, n * siglen * sizeof(int));

	cudaMemcpy(drowptr, &rowptr[0], rowptr.size() * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dcolidx, &colidx[0], colidx.size() * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dperms, hperms, n * siglen * sizeof(int), cudaMemcpyHostToDevice);

	getSig<<< (n + 1023 / 1024), 1024>>> (drowptr, dcolidx, dperms, dsigs, siglen, n);

	int sigs[n * siglen];
	cudaMemcpy(sigs, dsigs, n * siglen * sizeof(int), cudaMemcpyDeviceToHost);
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

	int num_bands = siglen / bsize;
	int *dbucks;
	cudaMalloc(&dbucks, n * num_bands * sizeof(int));

	getBuckets<<<(n + 1023) / 1024, 1024>>>(dsigs, dbucks, n, siglen, bsize, numbuckets);
	int hbucks[n * num_bands];
	cudaMemcpy(hbucks, dbucks, n * num_bands * sizeof(int), cudaMemcpyDeviceToHost);

	vector<set<int>> buckets(numbuckets);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < num_bands; j++) {
			int idx = hbucks[i * num_bands + j];
			buckets[idx].insert(i);
		}
	}

	set<pairi > result;
	for (auto s: buckets) {
		vi temp(s.begin(), s.end());
		for (int i = 0; i < temp.size(); i++) {
			for (int j = i + 1; j < temp.size(); j++) {
				result.insert(make_pair(temp[i], temp[j]));
			}
		}
	}
	return result;
}


class trio {
 public:
	int first, second, third;

	void print() {
		cout << first << " " << second << " " << third << endl;
	}

	trio(){}

	trio(int a, int b, int c) {
		first = a;
		second = b;
		third = c;
	}

};
class compare {
public:
   bool operator()(const trio &a, const trio &b) const {
	   if (a.f == b.f){
		   if(a.s==b.s)
			   return a.t<b.t;
		   return a.s<b.s;
	   }
	   return a.f > b.f;
   }
};

float J(vi &rowptr, vi &colidx, int i, int j){
	float ans = 0.0;
	int count=0;

	set<int>s;

	for(int k=rowptr[i];k<rowptr[i+1];++k){
		s.insert(colidx[k]);
		count++;
	}

	for(int k=rowptr[j];k<rowptr[j+1];++k){
		if(s.find(colidx[k])!=s.end()){
			ans++;
		}
		else{
			count++;
		}
	}


	return ans/count;
}

class mkDSU{
 public:
	vector<int>id,size,deleted;
	int nclusters,threshold_size;

	mkDSU(int n,int threshold){
		id.resize(n);
		size.resize(n);
		deleted.resize(n);

		for(int i=0; i<n; i++){
			id[i]=i;
			size[i]=1;
			deleted[i]=0;
		}
		nclusters=n-1;
		threshold_size=threshold;
	}

	int find(int a){
		int p=a,t;
		while(id[p]!=p)
			p=id[p];

		while(p!=a){
			t=id[a];
			id[a]=p;
			a=t;
		}
		return a;
	}

	void union_(set<pairi>& candidate_pairs, vi &rowptr, vi &colidx){

		set<trio,compare()> sim_queue;

		for(auto it:candidate_pairs)
			sim_queue.insert({J(rowptr,colidx,it.f,it.s),it.f,it.s});

		while(sim_queue.size() && nclusters>0){
			auto it = sim_queue.begin();
			trio temp =  *it;
			sim_queue.erase(it);

			int i = temp.s;
			int j = temp.t;

			if(i==id[i] && j==id[j]){
				if(deleted[i] || deleted[j])
					continue;

				if(size[i]<size[j]){
					id[i] = j;
					nclusters--;

					size[j] += size[i];

					if(size[j]>=threshold_size){
						deleted[j]=true;
						nclusters--;
					}
				}
				else{
					id[j] = i;
					clusters--;

					size[i] += size[j];

					if(size[i]>=threshold_size){
						deleted[i] = true;
						nclusters--;
					}
				}
			}
			else{ 

				int c1 = find(temp.s);
				int c2 = find(temp.t);

				if(deleted[c1] || deleted[c2] || c1==c2)
					continue;

				if(candidate_pairs.find({temp.s,temp.t})==candidate_pairs.end()){
					sim_queue.insert({J(S,c1,c2),min(c1,c2),max(c1,c2)});
					candidate_pairs.insert({min(c1,c2),max(c1,c2)});
				}
			}
		}
	}

	vi order_clusters(){
		map<int,vi> clusters;

		for(int i=0;i<n;++i){
			clusters[find(i)].push_back(i);
		}

		vi ans;

		for(auto it:clusters){
			for(auto ut:it){
				ans.push_back(ut);
			}
		}

		return ans;
	}
};





vi reordered_rows(vi &rowptr, vi &colidx){
	

	set<pairi> candidate_pairs = LSH(rowptr,colidx,5,5,5);

	mkDSU dsu(n,10);

	dsu.union_(candidate_pairs,rowptr,colidx);

	vi ans = dsu.order_clusters();

	return ans;
}






__global__ void MM(int *row_ptr, int *col_idx, int *col_val, int *dm, int *O, int N, int M, int K) {
	int idx = threadIdx.x + blockDim.x * blockIdx.x;

	if (idx < N * K) {
		// i'th row, j'th column of output
		int i = idx / K;
		int j = idx % K;
		int res = 0;
		int temp;
		for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
			temp = dm[col_idx[k] * K + j];
			res += col_val[k] * temp;
		}
		O[i * K + j] = res;
		// printf("%d %d %d\n", i, j, res);
	}
}

__global__ void SPMM(int *tile_row_ptr, int *panel_ptr, int *col_idx, int *col_val, int *D, int *O) {

	int row_panel_id = blockIdx.x;
	int row_id = threadIdx.x / 32;
	int thread_no = threadIdx.x % 32;

	int num_tiles = panel_ptr[row_panel_id + 1] - panel_ptr[row_panel_id];
	int global_row = PANEL_SIZE * row_panel_id + row_id;
	int ptr = panel_ptr[row_panel_id] * PANEL_SIZE + row_id * num_tiles;
	// printf("%d %d %d %d %d\n", row_panel_id, row_id, thread_no, num_tiles, ptr);

	for (int i = 0; i < num_tiles; ++i) {

		int low = tile_row_ptr[ptr + i];
		int high = tile_row_ptr[ptr + i + 1];

		for (int j = low; j < high; j++) {
			int temp = D[col_idx[j] * 32 + thread_no];
			O[global_row * 32 + thread_no] += col_val[j] * temp;
		}
	}
}

__global__ void find_dense_GPU(int *col_ptr, int *row_idx, int *isdense, int nr, int nc) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < nc) {
		int panel_id;
		for (int i = col_ptr[idx]; i < col_ptr[idx + 1]; ++i) {
			panel_id = row_idx[i] / PANEL_SIZE;
			isdense[panel_id * nc + idx] += 1;
		}
		// __syncthreads();
		// printf("%d %d %d\n", idx, counter, panel_id);
		int num_panels = nr / PANEL_SIZE;
		for (int i = 0; i < num_panels; i++) {
			int id = i * nc + idx;
			// printf("%d %d %d\n", idx, i, isdense[id]);
			if (isdense[id] >= DENSE_THRESHOLD)
				isdense[id] = 1;
			else
				isdense[id] = 0;
		}
	}
}

ve<vi > find_dense_CPU(vi &col_ptr, vi &row_idx, int nr, int nc) {
	int num_panels = nr / PANEL_SIZE;
	ve<vi > result(num_panels, vi());

	for (int j = 0; j < nc; j++) {
		int panel_id;
		vi is_dense(num_panels);
		for (int i = col_ptr[j]; i < col_ptr[j + 1]; i++) {
			panel_id = row_idx[i] / PANEL_SIZE;
			is_dense[panel_id]++;
		}
		for (int i = 0; i < num_panels; i++) {
			if (is_dense[i] >= DENSE_THRESHOLD) {
				result[i].push_back(j);
			}
		}
	}
	return result;
}

__global__ void ASPT_dense(int *tile_row_ptr,
									int *panel_ptr,
									int *col_idx,
									int *col_val,
									int *col_map,
									int *D,
									int *O) {

	__shared__ int shared_D[8192];

	int row_panel_id = blockIdx.x;
	int row_id = threadIdx.x / 32;
	int thread_no = threadIdx.x % 32;

	int num_tiles = panel_ptr[row_panel_id + 1] - panel_ptr[row_panel_id];
	int global_row = row_panel_id * PANEL_SIZE + row_id;
	int ptr = panel_ptr[row_panel_id] * PANEL_SIZE + row_id * num_tiles;

	for (int i = 0; i < num_tiles - 1; ++i) {

		int low = tile_row_ptr[ptr + i];
		int high = tile_row_ptr[ptr + i + 1];

		if (high > low) {
			shared_D[col_map[low] * 32 + thread_no] = D[col_idx[low] * 32 + thread_no];
		}

	}

	__syncthreads();

	for (int i = 0; i < num_tiles - 1; ++i) {

		int low = tile_row_ptr[i + ptr];
		int high = tile_row_ptr[i + ptr + 1];

		if (high > low) {

			int ind = col_map[low];
			O[global_row * 32 + thread_no] += col_val[low] * shared_D[ind * 32 + thread_no];
		}
	}

}

__global__ void ASPT_sparse(int *tile_row_ptr, int *panel_ptr, int *col_idx, int *col_val, int *D, int *O) {

	int row_panel_id = blockIdx.x;
	int row_id = threadIdx.x / 32;
	int thread_no = threadIdx.x % 32;

	int num_tiles = panel_ptr[row_panel_id + 1] - panel_ptr[row_panel_id];
	int global_row = row_panel_id * PANEL_SIZE + row_id;
	int ptr = panel_ptr[row_panel_id] * PANEL_SIZE + row_id * num_tiles + num_tiles - 1;

	int low = tile_row_ptr[ptr];
	int high = tile_row_ptr[ptr + 1];

	for (int i = low; i < high; ++i) {

		int j = col_idx[i];

		O[global_row * 32 + thread_no] += col_val[j] * D[j * 32 + thread_no];
	}

}

void run_MM(vi &row_ptr, vi &col_idx, vi &col_val, vi &host_DM, int nr, int nc, int ne) {
	// try a simple MM, (NxM)*(MxK) -> (NxK), use N*K threads

	int *drow_ptr;
	int *dcol_idx;
	int *dcol_val;
	cudaMalloc(&drow_ptr, (nr + 1) * sizeof(int));
	cudaMalloc(&dcol_idx, ne * sizeof(int));
	cudaMalloc(&dcol_val, ne * sizeof(int));

	cudaMemcpy(drow_ptr, &row_ptr[0], (nr + 1) * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dcol_idx, &col_idx[0], ne * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dcol_val, &col_val[0], ne * sizeof(int), cudaMemcpyHostToDevice);

	int *DM;
	cudaMalloc(&DM, nc * 32 * sizeof(int));
	cudaMemcpy(DM, &host_DM[0], nc * 32 * sizeof(int), cudaMemcpyHostToDevice);

	int *O;
	cudaMalloc(&O, nr * 32 * sizeof(int));
	MM <<< (nr * 32 + 1023) / 1024, 1024>>>(drow_ptr, dcol_idx, dcol_val, DM, O, nr, nc, 32);

	vi host_O(nr * 32);
	cudaMemcpy(&host_O[0], O, nr * 32 * sizeof(int), cudaMemcpyDeviceToHost);

	cudaFree(drow_ptr);
	cudaFree(dcol_idx);
	cudaFree(dcol_val);
	cudaFree(O);

	// for(int i=0; i<nr; i++)
	// {
	// 	for(int j=0; j<32; j++)
	// 		cout << host_O[32*i + j] << " ";
	// 	cout << endl;
	// }
	// cout << endl;
}

void run_SPMM(vi &tile_row_ptr, vi &panel_ptr, vi &col_idx, vi &col_val, vi &host_DM, int nr, int nc, int ne) {
	// trying SPMM with tiling (no reordering)
	int num_panels = nr / PANEL_SIZE;

	int *dtile_row_ptr;
	int *dpanel_ptr;
	int *dcol_idx;
	int *dcol_val;
	cudaMalloc(&dtile_row_ptr, tile_row_ptr.size() * sizeof(int));
	cudaMalloc(&dpanel_ptr, panel_ptr.size() * sizeof(int));
	cudaMalloc(&dcol_idx, ne * sizeof(int));
	cudaMalloc(&dcol_val, ne * sizeof(int));

	cudaMemcpy(dtile_row_ptr, &tile_row_ptr[0], tile_row_ptr.size() * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dpanel_ptr, &panel_ptr[0], panel_ptr.size() * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dcol_idx, &col_idx[0], ne * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dcol_val, &col_val[0], ne * sizeof(int), cudaMemcpyHostToDevice);
	int *DM;
	cudaMalloc(&DM, nc * 32 * sizeof(int));
	cudaMemcpy(DM, &host_DM[0], nc * 32 * sizeof(int), cudaMemcpyHostToDevice);

	int *O;
	cudaMalloc(&O, nr * 32 * sizeof(int));
	cudaMemset(O, 0, nr * 32 * sizeof(int));
	SPMM<<< num_panels, 32 * PANEL_SIZE>>>(dtile_row_ptr, dpanel_ptr, dcol_idx, dcol_val, DM, O);

	vi host_O(nr * 32);
	cudaMemcpy(&host_O[0], O, nr * 32 * sizeof(int), cudaMemcpyDeviceToHost);

	cudaFree(dtile_row_ptr);
	cudaFree(dpanel_ptr);
	cudaFree(dcol_idx);
	cudaFree(dcol_val);
	cudaFree(O);
	// for(int i=0; i<nr; i++)
	// {
	// 	for(int j=0; j<32; j++)
	// 		cout << host_O[32*i + j] << " ";
	// 	cout << endl;
	// }
	// cout << endl;
}

void run_ASPT(vi &tile_row_ptr,
				  vi &panel_ptr,
				  vi &col_idx,
				  vi &col_val,
				  vi &col_map,
				  vi &host_DM,
				  int nr,
				  int nc,
				  int ne) {
	// call ASPT kernels
	int num_panels = nr / PANEL_SIZE;

	int *dtile_row_ptr;
	int *dpanel_ptr;
	int *dcol_idx;
	int *dcol_val;
	int *dcol_map;
	cudaMalloc(&dtile_row_ptr, tile_row_ptr.size() * sizeof(int));
	cudaMalloc(&dpanel_ptr, panel_ptr.size() * sizeof(int));
	cudaMalloc(&dcol_idx, ne * sizeof(int));
	cudaMalloc(&dcol_val, ne * sizeof(int));
	cudaMalloc(&dcol_map, ne * sizeof(int));

	cudaMemcpy(dtile_row_ptr, &tile_row_ptr[0], tile_row_ptr.size() * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dpanel_ptr, &panel_ptr[0], panel_ptr.size() * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dcol_idx, &col_idx[0], ne * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dcol_val, &col_val[0], ne * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dcol_map, &col_map[0], ne * sizeof(int), cudaMemcpyHostToDevice);

	int *DM;
	cudaMalloc(&DM, nc * 32 * sizeof(int));
	cudaMemcpy(DM, &host_DM[0], nc * 32 * sizeof(int), cudaMemcpyHostToDevice);

	int *O1;
	cudaMalloc(&O1, nr * 32 * sizeof(int));
	cudaMemset(O1, 0, nr * 32 * sizeof(int));
	int *O2;
	cudaMalloc(&O2, nr * 32 * sizeof(int));
	cudaMemset(O2, 0, nr * 32 * sizeof(int));
	// cudaMalloc(&O2, nr*32*sizeof(int));
	// cudaMemset(O2, 0, nr*32*sizeof(int));
	// cudaDeviceSetCacheConfig(ASPT_dense, cudaFuncCachePreferShared);

	cudaStream_t s1, s2;
	cudaStreamCreate(&s1);
	cudaStreamCreate(&s2);

	ASPT_dense<<< num_panels, PANEL_SIZE * 32, 0, s1>>>(dtile_row_ptr, dpanel_ptr, dcol_idx, dcol_val, dcol_map, DM, O1);
	ASPT_sp
	arse<<<num_panels, PANEL_SIZE * 32, 0, s2>>>(dtile_row_ptr, dpanel_ptr, dcol_idx, dcol_val, DM, O2);

	vi host_O1(nr * 32);
	vi host_O2(nr * 32);
	cudaMemcpy(&host_O1[0], O1, nr * 32 * sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(&host_O2[0], O2, nr * 32 * sizeof(int), cudaMemcpyDeviceToHost);

	cudaFree(dtile_row_ptr);
	cudaFree(dpanel_ptr);
	cudaFree(dcol_idx);
	cudaFree(dcol_val);
	cudaFree(dcol_map);
	cudaFree(O1);
	cudaFree(O2);
	for (int i = 0; i < nr; i++) {
		for (int j = 0; j < 32; j++)
			cout << host_O1[32 * i + j] + host_O2[32 * i + j] << " ";
		cout << endl;
	}
}

void CSR_reorder_GPU(vi &rows, vi &cols, vi &row_ptr, vi &col_idx, vi &col_val, int nr, int nc, int ne) {
	// // create column wise CSR
	thrust::sort_by_key(cols.begin(), cols.begin()+ne, rows.begin());
	cout << "sorted cols" << endl;
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

	// find dense tiles now
	int num_panels = nr/PANEL_SIZE;
	int thr = num_panels*nc;
	cout << "tiles - " << thr << endl;

	int *dcol_ptr;
	int *drow_idx;
	int *is_dense;
	cudaMalloc(&dcol_ptr, (nc+1)*sizeof(int));
	cudaMalloc(&drow_idx, ne*sizeof(int));
	cudaMalloc(&is_dense, thr*sizeof(int));

	cudaMemcpy(dcol_ptr, &col_ptr[0], (nc+1)*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(drow_idx, &row_idx[0], ne*sizeof(int), cudaMemcpyHostToDevice);

	find_dense_GPU <<<(nc + 1023)/1024, 1024>>>(dcol_ptr, drow_idx, is_dense, nr, nc);

	vi isdense(thr);
	cudaMemcpy(&isdense[0], is_dense, thr*sizeof(int), cudaMemcpyDeviceToHost);
	cudaFree(dcol_ptr);
	cudaFree(drow_idx);
	cudaFree(is_dense);

	int total = 0;
	for(int i=0; i<num_panels; i++)
	{
		for(int j=0; j<nc; j++)
		{
			// cout << isdense[i*nc + j] << " ";
			if(isdense[i*nc+j])
				total ++;
		}
		// cout << endl;
	}
	cout << "dense tiles - " << total << endl;
	// Create row wise CSR
	// same as in cpu
}

int main(int argc, char **argv) {
	if (argc < 2)
		return 1;
	char *inputfilename = argv[1];
	FILE *fp;
	fp = fopen(inputfilename, "r");

	int nr, nc, ne;
	fscanf(fp, "%d %d %d", &nr, &nc, &ne);

	int r, c;
	vi rows(ne, 0);
	vi cols(ne, 0);
	for (int i = 0; i < ne; i++) {
		fscanf(fp, "%d %d", &r, &c);

		rows[i] = r;
		cols[i] = c;
	}

	thrust::sort_by_key(rows.begin(), rows.begin() + ne, cols.begin());

	vi temp_row_ptr(nr + 1, 0);
	vi temp_col_idx(ne, 0);

	for (int i = 0; i < ne; i++) {
		// cout << rows[i] << " " << cols[i] << endl;
		temp_row_ptr[rows[i]]++;
		temp_col_idx[i] = cols[i] - 1;
	}

	for (int i = 0; i < nr; i++)
		temp_row_ptr[i + 1] += temp_row_ptr[i];


	vi order_rows = reordered_rows(temp_row_ptr,temp_col_idx);


	for(int i=0;i<ne;++i){
		rows[i] = order_rows[rows[i]-1]; 
		cout << rows[i] << " ";
	}
	cout << endl;

	
	return 0;

	thrust::sort_by_key(cols.begin(), cols.begin() + ne, rows.begin());
	cout << "sorted cols" << endl;
	vi col_ptr(nc + 1, 0);
	vi row_idx(ne, 0);
	for (int i = 0; i < ne; i++) {
		// cout << cols[i] << " " << rows[i] << endl;
		col_ptr[cols[i]]++;
		row_idx[i] = rows[i] - 1;
	}

	for (int i = 0; i < nc; i++)
		col_ptr[i + 1] += col_ptr[i];

	// for(int i=0; i<=nc; i++)
	// 	cout << col_ptr[i] << " ";
	// cout <<endl;

	// for(int i=0; i<ne; i++)
	// 	cout << row_idx[i] << " ";
	// cout << endl;

	// find dense tiles now
	int num_panels = nr / PANEL_SIZE;

	ve<vi > dense = find_dense_CPU(col_ptr, row_idx, nr, nc);
	// for(int i=0; i<num_panels; i++)
	// {
	// 	for(int j=0; j<dense[i].size(); j++)
	// 	{
	// 		cout << dense[i][j] << " ";
	// 	}
	// 	cout << endl;
	// }
	thrust::sort_by_key(rows.begin(), rows.begin() + ne, cols.begin());
	cout << "sorted row wise" << endl;

	vi row_ptr(nr + 1, 0);
	vi col_idx(ne, 0);
	vi col_val(ne, 1);
	vi col_map(ne, 0);

	for (int i = 0; i < ne; i++) {
		// cout << rows[i] << " " << cols[i] << endl;
		row_ptr[rows[i]]++;
		col_idx[i] = cols[i] - 1;
	}

	for (int i = 0; i < nr; i++)
		row_ptr[i + 1] += row_ptr[i];

	cout << "Reordering" << endl;
	vi panel_ptr(num_panels + 1, 0);
	vi tile_row_ptr(1, 0);

	for (int panel_id = 0; panel_id < num_panels; ++panel_id) {
		map<int, int> densecols;
		for (int j = 0; j < dense[panel_id].size(); j++) {
			densecols[dense[panel_id][j]] = j;
		}

		panel_ptr[panel_id + 1] = densecols.size() + 1; // one sparse panel

		for (int i = panel_id * PANEL_SIZE; i < (panel_id + 1) * PANEL_SIZE; ++i) {
			if (i >= nr)
				break;

			ve<pairi > temp1;
			ve<pairi > temp2;
			for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
				if (densecols.find(col_idx[k]) == densecols.end()) {
					temp2.push_back(make_pair(col_idx[k], col_val[k]));
				} else {
					temp1.push_back(make_pair(col_idx[k], col_val[k]));
				}
			}

			int counter = 0;
			for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
				if (counter < temp1.size()) // dense columns
				{
					col_idx[k] = temp1[counter].f;
					col_val[k] = temp1[counter].s;
					col_map[k] = densecols[col_idx[k]];
				} else   // sparse columns
				{
					col_idx[k] = temp2[counter - temp1.size()].f;
					col_val[k] = temp2[counter - temp1.size()].s;
					col_map[k] = -1;
				}
				++counter;
			}

			counter = 0;
			int found = 0;
			for (auto mapel:densecols) {
				auto el = mapel.f;
				// cout << el << ' ';
				found = 0;
				for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
					if (el == col_idx[k]) {
						found = 1;
						counter++;
						break;
					} else if (el < col_idx[k]) {
						break;
					}
				}

				tile_row_ptr.push_back(found);

				// cout << el << " " << tile_row_ptr[tile_row_ptr.size()-1] << found << endl;
			}
			// cout << endl;
			tile_row_ptr.push_back(row_ptr[i + 1] - row_ptr[i] - counter);

		}

		// densecols.clear();
	}

	for (int i = 0; i < num_panels; i++)
		panel_ptr[i + 1] += panel_ptr[i];

	for (int i = 1; i < tile_row_ptr.size(); i++)
		tile_row_ptr[i] += tile_row_ptr[i - 1];

	// cout << "row_ptr" << endl;
	// for(int i=0; i<=nr; i++)
	// 	cout << row_ptr[i] << " ";
	// cout <<endl;

	// cout << "col_idx" << endl;
	// for(int i=0; i<ne; i++)
	// 	cout << col_idx[i] << " ";
	// cout << endl;

	// cout << "panel_ptr" << endl;
	// for(int i=0; i<= num_panels; i++)
	// 	cout << panel_ptr[i] << " ";
	// cout << endl;
	// cout << "tile_row_ptr" << endl;
	// for(int i=0; i<tile_row_ptr.size(); i++)
	// 	cout << tile_row_ptr[i] << " ";
	// cout << endl;
	// cout << endl;

	// cout << "dense col mapping" << endl;
	// for(int i=0; i<col_map.size(); i++)
	// 	cout << col_map[i] << " ";
	// cout << endl;


	vi host_DM(nc * 32, 1);

	int *DM;
	cudaMalloc(&DM, nc * 32 * sizeof(int));
	cudaMemcpy(DM, &host_DM[0], nc * 32 * sizeof(int), cudaMemcpyHostToDevice);

	// run_MM(row_ptr, col_idx, col_val, host_DM, nr, nc, ne);
	// run_SPMM(tile_row_ptr, panel_ptr, col_idx, col_val, host_DM, nr, nc, ne);
	// run_ASPT(tile_row_ptr, panel_ptr, col_idx, col_val, col_map, host_DM, nr, nc, ne);

}

