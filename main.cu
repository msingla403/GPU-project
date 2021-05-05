#include<bits/stdc++.h>
using namespace std;

#define pairi pair<int,int>
#define ve vector
#define vi vector<int>
#define f first
#define s second
#define t third

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

};

float J(ve<vi>&S, int i, int j){
	float ans = 0.0;
	int count=0;

	for(int k=0;k<S[i].size();++k){
		if(S[i][k]!=0 && S[j][k]!=0){
			ans++;
		}
		if(S[i][k]!=0 || S[j][k]!=0){
			count++;
		}
	}
	return ans/count;
}

set<pairi> LSH(ve<vi>&S, int siglen, int bsize){


}




class mkDSU{
 public:
	vector<int>id,size,deleted,nclusters,threshold_size;

	mkDSU(int n,int threshold){
		id.resize(n);
		size.resize(n);
		deleted.resize(n);

		FOR(i,0,n){
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

	void union_(set<pairi>& candidate_pairs, ve<vi>& S){

		set<trio,trio::compare()> sim_queue;

		for(auto it:candidate_pairs)
			sim_queue.insert({J(S,it.f,it.s),it.f,it.s});

		while(sim_queue.size() && nclusters>0){
			auto it = sim_queue.begin();
			trio temp =  *it;
			sim_queue.erase(it);

			
			

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





vi reordered_rows(ve<vi>&S){
	int n = S.size();

	set<pairi> candidate_pairs = LSH(S,5,5);


	mkDSU dsu(n);

	dsu.union_(candidate_pairs);

	vi ans = dsu.order_clusters();

	return ans;
}


#define PANEL_SIZE 3

__global__ void SPMM(int * tile_row_ptr, int * panel_ptr, int * col_val, int * col_idx){

	int row_panel_id = blockIdx.x;
	int row_id = threadIdx.x/32;
	int thread_no = threadIdx.x%32;

	int num_tiles = panel_ptr[row_panel_id+1] - panel_ptr[row_panel_id];

	int ptr = panel_ptr[row_panel_id]*PANEL_SIZE + row_id*num_tiles;

	for(int i=0;i<num_tiles;++i){

		int low = tile_row_ptr[ptr+i];
		int high = tile_row_ptr[ptr+i+1];

		if(high>low){
			int j=low;
			O[row_id][thread_no] += col_val[j] * D[col_idx[j]][thread_no];
		}
	}
}

__global__ void ASPT_dense(int * panel_ptr, int * col_val, int * col_idx ){

	int row_panel_id = blockIdx.x;
	int row_id = threadIdx.x/32;
	int thread_no = threadIdx.x%32;

	int num_tiles = panel_ptr[row_panel_id+1] - panel_ptr[row_panel_id];

	int ptr = panel_ptr[row_panel_id]*PANEL_SIZE + row_id*num_tiles;

	__shared__ int map_tiles[(num_tiles-1)*PANEL_SIZE];
	__shared__ int shared_D[num_tiles-1][32];

	if(thread_no==0){
		for(int i=0;i<num_tiles-1;++i){

			int low = tile_row_ptr[ptr+i];
			int high = tile_row_ptr[ptr+i+1];

			if(high>low){
				map_tiles[i]=col_idx[low];
			}
		}

	}

	__syncthreads();


	

	__syncthreads();

	for(int i=0;i<num_tiles;++i){

		int low = tile_row_ptr[i+ptr];
		int high = tile_row_ptr[i+ptr+1];


		for(int j=low;j<=high;++j){
			O[row_id][thread_no] += col_val[j] * D[col_idx[j]][thread_no];
		}
	}
}









