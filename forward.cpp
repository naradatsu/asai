#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

using namespace std;

void convert_to_num(string s, int* S, int a, int alphabet_size, string alphabet){
    for(int i = 0; i < a; i++){
        for(int j = 0; j < alphabet_size; j++){
            if(s[i] == alphabet[j]){
                S[i] = j;
                break;
            }
        }
    }
}

int main(void){
    ifstream ist("hmminfo.txt");
    if(!ist){
        cout << "can't open" << endl;
        exit(1);
    }

    int alphabet_size = 0;
    string alphabet;
    int state_num = 0;
    double* state_transision_prob; //1 dimension array of state_num * state_num. a_ij is given by [i*state_num + j]
    double* output_prob; //1 dimension array of state_num * aiphabet_size. e(i, s) is given by [(i-1)*state_num + (order of s)]

    ist >> alphabet_size;

    for(int i = 0; i < alphabet_size; i++){
        string p;
        ist >> p;
        alphabet += p;
    }

    ist >> state_num;
    state_transision_prob = (double*)malloc(sizeof(double) * (state_num * state_num));
    output_prob = (double*)malloc(sizeof(double) * ((state_num - 1) * alphabet_size));

    int a = state_num * state_num;
    int b = (state_num - 1) * alphabet_size;

    for(int i = 0; i < a; i++){
        ist >> state_transision_prob[i];
        state_transision_prob[i] = log(state_transision_prob[i]);
        //cout << state_transision_prob[i] << " ";
    }

    //cout << "output" <<endl;
    for(int i = 0; i < b; i++){
        ist >> output_prob[i];
        output_prob[i] = log(output_prob[i]);
        //cout << output_prob[i] << " ";
    }

    ist.close();

    ifstream isf("observed_seq.fasta");
    if(!isf){
        cout << "can't open" << endl;
        exit(1);
    }

    string s;
    string c; 
    getline(isf, c);
    while(!isf.eof()){
        string p;
        getline(isf, p);
        s += p;
    }

    //cout << s;
    
    int arraysize = s.size();
    int* S = (int*)malloc(arraysize * sizeof(int));
    convert_to_num(s, S, arraysize, alphabet_size, alphabet);

    //cout << arraysize;

    
    //for(int i = 0; i < arraysize; i++) cout << S[i] << " ";

    double* DP_table = (double*)malloc(sizeof(double) * (state_num - 1) * arraysize); //in order to access v_ij (i:state, j:letter), [(j)*(state_num) + i].
    int* origin_table = (int*)malloc(sizeof(int) * (state_num - 1) * arraysize); // for traceback.

    //cout << endl;
    //initialize
    for(int i = 1; i < state_num; i++){
        //cout << "stp[" << i <<"] + opp[" << (i-1)*alphabet_size + S[0] << "]  ";
        DP_table[i - 1] = state_transision_prob[i] + output_prob[(i-1)*alphabet_size + S[0]];
        origin_table[i - 1] = 0;
        //cout << DP_table[i - 1] << endl;
    }

    //cout << endl;

    //<arraysize
    for(int i = 1; i < arraysize; i++){
        for(int j = 1; j < state_num; j++){
            //find max
            double* candidates = (double*)malloc(sizeof(double) * state_num);
            for(int k = 1; k < state_num; k++){
                //cout << "DPtable["<<(i - 1) * (state_num - 1) + (k - 1) << "], " << "stprob[" << k * state_num + j << "], otprob[" << (j-1) * alphabet_size + S[i]<< "]    ";
                candidates[k] = DP_table[(i - 1) * (state_num - 1) + (k - 1)] + state_transision_prob[k * state_num + j] + output_prob[(j-1) * alphabet_size + S[i]];
            }
            int origin = 0;
            double biggest = -100000000;
            for(int k = 1; k < state_num; k++){
                if(candidates[k] >= biggest){
                    biggest = candidates[k];
                    origin = k;
                }
            }
            DP_table[i * (state_num - 1) + j - 1] = biggest;
            origin_table[i * (state_num - 1) + j - 1] = origin;

            //cout << "(" << i << "," << j << ")";
            //cout << DP_table[i * (state_num - 1) + j - 1] << "  ";
            //cout << "ori:" << origin;
            //cout << endl;
        }
    }

    int final_state = 0;
    double probability = -1000000000;
    for(int i = 1; i < state_num; i++){
        if(DP_table[(state_num - 1) * (arraysize - 1) + i - 1] >= probability){
            final_state = origin_table[(state_num - 1) * (arraysize - 1) + i - 1];
            probability = DP_table[(state_num - 1) * (arraysize - 1) + i - 1];
        }
    }

    int trace = final_state;
    int* route = (int*)malloc(sizeof(int) * arraysize);
    for(int i = arraysize - 1;i >= 0; i--){
        if(trace == 0) break;
        route[i] = trace;
        trace = origin_table[(state_num - 1) * i + trace - 1];
    }

    cout << "0";
    for(int i = 0; i < arraysize; i++){
        cout << route[i];
    }

    return 0;
}