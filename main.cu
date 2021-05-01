#include<bits/stdc++.h>
using namespace std;

#define pairi pair<int,int>
#define ve vector
#define vi vector<int>


class mkDSU{
 public:
	vector<int>par,size;

	mkDSU(int n){
		par.resize(n);
		size.resize(n);

		FOR(i,0,n){
			par[i]=i;
			size[i]=0;
		}
	}

	int find(int a){
		int p=a,t;
		while(par[p]!=p)
			p=par[p];

		while(p!=a){
			t=par[a];
			par[a]=p;
			a=t;
		}
		return a;
	}

	void union_(int a,int b){
		int u1=find(a);
		int u2=find(b);

		if(u1 == u2)
			return ;

		if(size[u1]<size[u2]){
			par[u1]=u2;
		}else{
			if(size[u2]<size[u1]){
				par[u2]=u1;
			}
			else{
				par[u1]=u2;
				size[u2]++;
			}
		}
	}
};


ve<pairi> LSH(ve<vi>&S, int siglen, int bsize){


}

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

vi reordered_rows(ve<vi>&S){
	int n = S.size();
}