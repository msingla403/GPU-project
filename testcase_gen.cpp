#include <iostream>
#include <algorithm>
#include <vector>
using namespace std;

int main()
{
    srand(time(NULL));
    long int nr = 1000;
    long int nc = 1000;
    float sparsity = 0.05;

    double ne = 1.0*sparsity*nr*nc;
    // cout << ne << endl;
    int count = 0;

    vector< pair<int,int> > posi;

    vector<int> perm(nr);
    for(int i=0; i<nr; i++)
        perm[i] = i;
    
    for(int i=0; i<nr; i++)
    {
        random_shuffle(perm.begin(), perm.end());

        int var = -20 + rand()%40;
        float el = 1.0*sparsity*(100+var)/100;
        el *= nc;
        // cout << el << endl;
        for(int j=0; j< el; j++)
        {
            posi.push_back(make_pair(i+1, perm[j]+1));
            count++;
            if(count >= ne)
                break;
        }
        if(count >= ne)
            break;
    }

    cout << nr << " " << nc << " " << count << endl;
    for(int i=0; i<posi.size(); i++)
    {
        cout << posi[i].first << " " << posi[i].second << endl;
    }
}