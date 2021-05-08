#include<bits/stdc++.h>
#include<cuda.h>
using namespace std;

#define pairi pair<int,int>
#define ve vector
#define vi vector<int>
#define f first
#define s second
#define t third

// class trio {
//  public:
// 	int first, second, third;

// 	void print() {
// 		cout << first << " " << second << " " << third << endl;
// 	}

// 	trio(){}

// 	trio(int a, int b, int c) {
// 		first = a;
// 		second = b;
// 		third = c;
// 	}

// 	class compare {
// 	 public:
// 		bool operator()(const trio &a, const trio &b) const {
// 			if (a.f == b.f){
// 				if(a.s==b.s)
// 					return a.t<b.t;
// 				return a.s<b.s;
// 			}
// 			return a.f > b.f;
// 		}
// 	};

// };

// float J(ve<vi>&S, int i, int j){
// 	float ans = 0.0;
// 	int count=0;

// 	for(int k=0;k<S[i].size();++k){
// 		if(S[i][k]!=0 && S[j][k]!=0){
// 			ans++;
// 		}
// 		if(S[i][k]!=0 || S[j][k]!=0){
// 			count++;
// 		}
// 	}
// 	return ans/count;
// }

// class mkDSU{
//  public:
// 	vector<int>id,size,deleted,nclusters,threshold_size;

// 	mkDSU(int n,int threshold){
// 		id.resize(n);
// 		size.resize(n);
// 		deleted.resize(n);

// 		FOR(i,0,n){
// 			id[i]=i;
// 			size[i]=1;
// 			deleted[i]=0;
// 		}
// 		nclusters=n-1;
// 		threshold_size=threshold;
// 	}

// 	int find(int a){
// 		int p=a,t;
// 		while(id[p]!=p)
// 			p=id[p];

// 		while(p!=a){
// 			t=id[a];
// 			id[a]=p;
// 			a=t;
// 		}
// 		return a;
// 	}

// 	void union_(set<pairi>& candidate_pairs, ve<vi>& S){

// 		set<trio,trio::compare()> sim_queue;

// 		for(auto it:candidate_pairs)
// 			sim_queue.insert({J(S,it.f,it.s),it.f,it.s});

// 		while(sim_queue.size() && nclusters>0){
// 			auto it = sim_queue.begin();
// 			trio temp =  *it;
// 			sim_queue.erase(it);

			
			

// 			if(i==id[i] && j==id[j]){
// 				if(deleted[i] || deleted[j])
// 					continue;

// 				if(size[i]<size[j]){
// 					id[i] = j;
// 					nclusters--;

// 					size[j] += size[i];

// 					if(size[j]>=threshold_size){
// 						deleted[j]=true;
// 						nclusters--;
// 					}

// 				}
// 				else{
// 					id[j] = i;
// 					clusters--;

// 					size[i] += size[j];

// 					if(size[i]>=threshold_size){
// 						deleted[i] = true;
// 						nclusters--;
// 					}

// 				}

// 			}


// 			else{ 

// 				int c1 = find(temp.s);
// 				int c2 = find(temp.t);

// 				if(deleted[c1] || deleted[c2] || c1==c2)
// 					continue;

// 				if(candidate_pairs.find({temp.s,temp.t})==candidate_pairs.end()){
// 					sim_queue.insert({J(S,c1,c2),min(c1,c2),max(c1,c2)});
// 					candidate_pairs.insert({min(c1,c2),max(c1,c2)});
// 				}
// 			}

// 		}
		
// 	}

// 	vi order_clusters(){
// 		map<int,vi> clusters;

// 		for(int i=0;i<n;++i){
// 			clusters[find(i)].push_back(i);
// 		}

// 		vi ans;

// 		for(auto it:clusters){
// 			for(auto ut:it){
// 				ans.push_back(ut);
// 			}
// 		}

// 		return ans;
// 	}
// };





// vi reordered_rows(ve<vi>&S){
// 	int n = S.size();

// 	set<pairi> candidate_pairs = LSH(S,5,5);


// 	mkDSU dsu(n);

// 	dsu.union_(candidate_pairs);

// 	vi ans = dsu.order_clusters();

// 	return ans;
// }


// #define PANEL_SIZE 3

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

// __global__ void ASPT_dense(int * panel_ptr, int * col_val, int * col_idx ){

// 	int row_panel_id = blockIdx.x;
// 	int row_id = threadIdx.x/32;
// 	int thread_no = threadIdx.x%32;

// 	int num_tiles = panel_ptr[row_panel_id+1] - panel_ptr[row_panel_id];

// 	int ptr = panel_ptr[row_panel_id]*PANEL_SIZE + row_id*num_tiles;

// 	__shared__ int map_tiles[(num_tiles-1)*PANEL_SIZE];
// 	__shared__ int shared_D[num_tiles-1][32];

// 	if(thread_no==0){
// 		for(int i=0;i<num_tiles-1;++i){

// 			int low = tile_row_ptr[ptr+i];
// 			int high = tile_row_ptr[ptr+i+1];

// 			if(high>low){
// 				map_tiles[i]=col_idx[low];
// 			}
// 		}

// 	}

// 	__syncthreads();


	

// 	__syncthreads();

// 	for(int i=0;i<num_tiles;++i){

// 		int low = tile_row_ptr[i+ptr];
// 		int high = tile_row_ptr[i+ptr+1];


// 		for(int j=low;j<=high;++j){
// 			O[row_id][thread_no] += col_val[j] * D[col_idx[j]][thread_no];
// 		}
// 	}
// }

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

int main()
{
	// int n = 6;
	// int m = 6;
	vi rowptr{0,2,5,7,8,11,13};
	vi colidx{0,4,1,3,5,2,4,1,0,3,4,2,5};
	
	set<pairi> candidates = LSH(rowptr, colidx, 6, 2, 6);

	for(auto i:candidates)
	{
		cout << i.f << " " << i.s << endl;
	}
}



