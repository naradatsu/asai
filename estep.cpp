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

//i:state, alphabet.
double acE(int i, int alphabet, double* E, int alphabet_size){
    return E[(i-1)*alphabet_size + alphabet];
}

//from state i to j.
double acA(int i, int j, double* A, int state_num){
    return A[i*state_num + j];
}

double* acpE(int i, int alphabet, double* E, int alphabet_size){
    return &E[(i-1)*alphabet_size + alphabet];
}

double* acpA(int i, int j, double* A, int state_num){
    return &A[i*state_num + j];
}

double acdp(int i, int t, double* DP, int state_num){
    return DP[(t-1)*(state_num-1) + (i-1)];
}

double* acpdp(int i, int t, double* DP, int state_num){
    return &DP[(t-1)*(state_num-1) + (i-1)];
}

void showdp(double* dp, int arraysize, int state_num){
    for(int i = 1; i < state_num; i++){
        cout << i << " : ";
        for(int t = 1; t < arraysize + 1; t++){
            cout << acdp(i, t, dp, state_num) << "   ";
        }
        cout << endl;
    }
}

int main(void){
    //0. PREPARE
    ifstream ist("testinfo.txt");
    if(!ist){
        cout << "can't open" << endl;
        exit(1);
    }

    int alphabet_size = 0;
    string alphabet;
    int state_num = 0;
    double* A; //1 dimension array of state_num * state_num. a_ij is given by [i*state_num + j]
    double* E; //1 dimension array of state_num * aiphabet_size. e(i, s) is given by [(i-1)*state_num + (order of s)]

    ist >> alphabet_size;

    for(int i = 0; i < alphabet_size; i++){
        string p;
        ist >> p;
        alphabet += p;
    }

    ist >> state_num;
    A = (double*)malloc(sizeof(double) * (state_num * state_num));
    E = (double*)malloc(sizeof(double) * ((state_num - 1) * alphabet_size));

    for(int i = 0; i < state_num * state_num; i++){
        ist >> A[i];
        A[i] = A[i];
    }

    for(int i = 0; i < (state_num - 1) * alphabet_size; i++){
        ist >> E[i];
        E[i] = E[i];
    }

    ist.close();

    ifstream isf("test.fasta");
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

    int arraysize = s.size();
    int* S = (int*)malloc(arraysize * sizeof(int));
    convert_to_num(s, S, arraysize, alphabet_size, alphabet);

    //1. FORWARD ALGORITHM
    double* dpfor = (double*)malloc(sizeof(double) * (state_num - 1) * arraysize);
    //initialize f_l(x1) 
    for(int l = 1; l < state_num; l++){
        dpfor[(1 - 1)*(state_num-1) + l - 1] = acE(l, S[0], E, alphabet_size) * acA(0, l, A, state_num);
    }

    //Recurrance
    for(int t = 2; t < arraysize+1; t++){
        for(int l = 1; l < state_num; l++){
            double sigma = 0;
            for(int k = 1; k < state_num; k++){
                sigma += acdp(k, t-1, dpfor, state_num) * acA(k, l, A, state_num);
            }
            //dpfor[(t - 1)*(state_num - 1) + l - 1] = acE(l, S[t-1], E, alphabet_size) * sigma;
            *acpdp(l, t, dpfor, state_num) = acE(l, S[t-1], E, alphabet_size) * sigma;
        }
    }

    //Finshing process
    double Px = 0;
    for(int l = 1; l < state_num; l++){
        Px += acdp(l, arraysize, dpfor, state_num);
    }

    //2. BACKWARD ALGORITHM
    double* dpback = (double*)malloc(sizeof(double) * (state_num - 1) * arraysize); 

    //initialize b_l(xT) 
    for(int l = 1; l < state_num; l++){
        dpback[(arraysize-1)*(state_num-1) + l - 1] = 1;
    }

    //Recurrence
    for(int t = arraysize - 2; t >= 0; t--){
        for(int k = 1; k < state_num; k++){
            double sigma = 0;
            for(int l = 1; l < state_num; l++){
                sigma += acA(k, l, A, state_num) * acE(l, S[t+1], E, alphabet_size) * acdp(l, t+2, dpback, state_num);
            }
            *acpdp(k, t+1, dpback, state_num) = sigma;
        }
    }

    //finishing process
    double Pxb = 0;
    for(int l = 1; l < state_num; l++){
        Pxb += acA(0, l, A, state_num) * acE(l, S[0], E, alphabet_size) * acdp(l, 1, dpback, state_num);
    }

    cout << "Px for : " << Px << endl << "Px back: " << Pxb << endl;

    

    //3. Caluculate merginalized prob. (M STEP)
    double* Ahat = (double*)malloc(sizeof(double) * (state_num * state_num));
    for(int i = 0; i < state_num; i++){
        for(int j = 0; j < state_num; j++){
            double sigma = 0;
            for(int t = 1; t < arraysize; t++){
                sigma += acdp(i, t, dpfor, state_num) * acA(i, j, A, state_num) * acE(j, S[t], E, alphabet_size) * acdp(j, t+1, dpback, state_num);
            }
            *acpA(i, j, Ahat, state_num) = sigma / Px;
        }
    }


    double* Ehat;
    Ehat = (double*)malloc(sizeof(double) * (state_num-1)*alphabet_size);
    for(int i = 1; i < state_num; i++){
        for(int c = 0; c < alphabet_size; c++){
            double sigma = 0;
            for(int t = 0; t < arraysize; t++){
                if(S[t] == c){
                    sigma += acdp(i, t+1, dpfor, state_num) * acdp(i, t+1, dpback, state_num);
                }
            }
            *acpE(i, c, Ehat, alphabet_size)= sigma / Px;
        }
    }

    double* Anew =  (double*)malloc(sizeof(double) * (state_num * state_num));
    double* Enew =  (double*)malloc(sizeof(double) * (state_num-1)*alphabet_size);
    for(int k = 0; k < state_num; k++){
        double sigma = 0;
        for(int i = 0; i < state_num; i++){
            sigma += acA(k, i, Ahat, state_num);
        }
        for(int l = 0; l < state_num; l++){
            *acpA(k, l, Anew, state_num) = acA(k, l, Ahat, state_num) / sigma;
        }
    }

    for(int k = 1; k < state_num; k++){
        double sigma = 0;
        for(int b = 0; b < alphabet_size; b++){
            sigma += acE(k, b, Ehat, alphabet_size);
        }
        for(int b = 0; b < alphabet_size; b++){
            *acpE(k, b, Enew, alphabet_size) = acE(k, b, Ehat, alphabet_size) / sigma;
        }
    }

    //1. FORWARD ALGORITHM
    double* dpfor2 = (double*)malloc(sizeof(double) * (state_num - 1) * arraysize);
    //initialize f_l(x1) 
    for(int l = 1; l < state_num; l++){
        dpfor2[(1 - 1)*(state_num-1) + l - 1] = acE(l, S[0], Enew, alphabet_size) * acA(0, l, Anew, state_num);
    }

    //Recurrance
    for(int t = 2; t < arraysize+1; t++){
        for(int l = 1; l < state_num; l++){
            double sigma = 0;
            for(int k = 1; k < state_num; k++){
                sigma += acdp(k, t-1, dpfor2, state_num) * acA(k, l, Anew, state_num);
            }
            //dpfor[(t - 1)*(state_num - 1) + l - 1] = acE(l, S[t-1], E, alphabet_size) * sigma;
            *acpdp(l, t, dpfor2, state_num) = acE(l, S[t-1], Enew, alphabet_size) * sigma;
        }
    }

    //Finshing process
    double Px2 = 0;
    for(int l = 1; l < state_num; l++){
        Px2 += acdp(l, arraysize, dpfor2, state_num);
    }

    cout << "new Px for : " << Px2;
    
    return 0;
}