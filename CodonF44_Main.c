// Author: Zhang Zetong
// Date: 2024-08-17
#define PY_SSIZE_T_CLEAN

#include<Python.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#define MAX_ID_LEN 256
#define MAX_LEN 2048
#define MAX_ITER 0xffff
#define MAX_COL 60

#define N -1
#define A  0
#define G  1
#define C  2
#define T  3

typedef struct CodonType {
    char seq[4];
    float cai;
    float freqs;
} Codon;

Codon CODONS[21][6] = {
/* C */	{{"TGT", 1.00, 0.79}, {"TGC", 0.26, 1.00}}, 
/* D */	{{"GAT", 1.00, 0.74}, {"GAC", 0.34, 1.00}}, 
/* S */	{{"TCA", 1.00, 0.34}, {"TCT", 0.76, 0.61}, {"AGT", 0.63, 0.82}, {"AGC", 0.23, 0.90}, {"TCG", 0.15, 0.95}, {"TCC", 0.13, 1.00}}, 
/* Q */	{{"CAA", 1.00, 0.84}, {"CAG", 0.19, 1.00}}, 
/* M */	{{"ATG", 1.00, 1.00}}, 
/* N */	{{"AAT", 1.00, 0.80}, {"AAC", 0.25, 1.00}}, 
/* P */	{{"CCA", 1.00, 0.48}, {"CCT", 0.76, 0.84}, {"CCG", 0.17, 0.92}, {"CCC", 0.16, 1.00}}, 
/* K */	{{"AAA", 1.00, 0.84}, {"AAG", 0.19, 1.00}}, 
/* * */	{{"TAA", 1.00, 0.73}, {"TGA", 0.22, 0.89}, {"TAG", 0.16, 1.00}}, 
/* T */	{{"ACA", 1.00, 0.40}, {"ACT", 0.93, 0.78}, {"ACC", 0.28, 0.89}, {"ACG", 0.27, 1.00}}, 
/* F */	{{"TTT", 1.00, 0.76}, {"TTC", 0.32, 1.00}}, 
/* A */	{{"GCT", 1.00, 0.42}, {"GCA", 0.77, 0.75}, {"GCC", 0.34, 0.90}, {"GCG", 0.24, 1.00}}, 
/* G */	{{"GGT", 1.00, 0.39}, {"GGA", 0.98, 0.77}, {"GGC", 0.31, 0.88}, {"GGG", 0.30, 1.00}}, 
/* I */	{{"ATT", 1.00, 0.70}, {"ATC", 0.27, 0.89}, {"ATA", 0.16, 1.00}}, 
/* L */	{{"TTA", 1.00, 0.33}, {"CTT", 0.78, 0.59}, {"TTG", 0.63, 0.80}, {"CTC", 0.22, 0.87}, {"CTA", 0.22, 0.95}, {"CTG", 0.15, 1.00}}, 
/* H */	{{"CAT", 1.00, 0.75}, {"CAC", 0.33, 1.00}}, 
/* R */	{{"CGT", 1.00, 0.42}, {"AGA", 0.54, 0.65}, {"CGA", 0.34, 0.80}, {"CGC", 0.26, 0.91}, {"CGG", 0.14, 0.97}, {"AGG", 0.09, 1.00}}, 
/* W */	{{"TGG", 1.00, 1.00}}, 
/* V */	{{"GTT", 1.00, 0.51}, {"GTA", 0.39, 0.71}, {"GTC", 0.32, 0.87}, {"GTG", 0.25, 1.00}}, 
/* E */	{{"GAA", 1.00, 0.83}, {"GAG", 0.20, 1.00}}, 
/* Y */	{{"TAT", 1.00, 0.78}, {"TAC", 0.27, 1.00}}, 
};

const int NUMS[21] = { 2, 2, 6, 1, 2, 2, 4, 2, 3, 4, 2, 4, 4, 3, 6, 2, 6, 1, 4, 2, 2 };

const char SITES[][7] = { "GAATTC", "ACTAGT", "TCTAGA", "GTCGAC" };

void coding(char *seq, double *params, int len) {
    double nt[4] = {0}, pt[3][4] = {{0}}, dt[4][4] = {{0}};
    int i, j, n;
    short arr[0xffff];

    for (i = 0; i < len; i++)
        switch(seq[i]) {
            case 'a':
            case 'A': arr[i] = A; break;
            case 'g':
            case 'G': arr[i] = G; break;
            case 'c':
            case 'C': arr[i] = C; break;
            case 't':
            case 'T': arr[i] = T; break;
            default : arr[i] = N;
        }
    
    for (i = 0, n = len; i < len; i++) {
        nt[arr[i]] ++;
        pt[i % 3][arr[i]] ++;
        if (i < len - 1) dt[arr[i]][arr[i + 1]] ++;
    }

    params[0] = (nt[A] + nt[G] - nt[C] - nt[T]) / len;
    params[1] = (nt[A] + nt[C] - nt[G] - nt[T]) / len;
    params[2] = (nt[A] + nt[T] - nt[C] - nt[G]) / len;
        
    for (i = 0, params += 3; i < 3; i++, params += 3) {
        params[0] = pt[i][A] + pt[i][G] - pt[i][C] - pt[i][T];
        params[1] = pt[i][A] + pt[i][C] - pt[i][G] - pt[i][T];
        params[2] = pt[i][A] + pt[i][T] - pt[i][C] - pt[i][G];
        
        n = pt[i][A] + pt[i][G] + pt[i][C] + pt[i][T];

        if (n > 0) for (j = 0; j < 3; j++) params[j] /= n;
    }

    for (i = 0; i < 4; i ++, params += 3) {
        params[0] = dt[i][A] + dt[i][G] - dt[i][C] - dt[i][T];
        params[1] = dt[i][A] + dt[i][C] - dt[i][G] - dt[i][T];
        params[2] = dt[i][A] + dt[i][T] - dt[i][C] - dt[i][G];

        n = pt[i][A] + pt[i][G] + pt[i][C] + pt[i][T];

        if (n > 0) for (j = 0; j < 3; j++) params[j] /= n;
    }
}

int main(int argc, char *argv[]) {

    FILE *in = fopen(argv[1], "r");

    if (in == NULL) {
        printf("[ERROR]Invalid input file '%s'. \n", argv[1]);
        exit(-1);
    }
    
    strcpy(strrchr(argv[1], '.'), ".out.fa");

    FILE *out = fopen(argv[1], "w");

    if (out == NULL) {
        printf("[ERROR]Fail to create output file '%s'. \n", argv[1]);
        exit(1);
    }

    char buf[MAX_ID_LEN], tri[4] = {0}, bst[MAX_LEN * 3], chr;
    int pro[MAX_LEN], aa, cdn, len, i, j, res;
    double ran, bfr, aft, cai, params[24];

    Py_Initialize();

    PyObject *SVM, *pModul, *vec = PyList_New(24);
    
    if (!Py_IsInitialized()) {
        printf("[ERROR]Fail to init Python. \n");
        exit(2);
    } else {
        PyRun_SimpleString("import sys\n");
        PyRun_SimpleString("sys.path.append('./')\n");
        
        pModul = PyImport_ImportModule("CodonF44_Load");

        if (pModul == NULL) {
            printf("[ERROR]Cannot find CodonF44_Load.py. \n");
            exit(3);
        }

        SVM = PyObject_GetAttrString(pModul, "SVM");

        if (SVM == NULL) {
            printf("[ERROR]Cannot find SVM Loader in CodonF44_Load.py. \n");
            exit(4);
        }
    }

    while (fgets(buf, sizeof(buf), in)) {
        i = len = bfr = 0;

        while (!feof(in) && (chr = fgetc(in)) != '>') {
            if (chr >= 'a' && chr <= 'z')
                chr -= 32;
            
            if (chr >= 'A' && chr <= 'Z')
                tri[i ++] = chr;
            
            if (i == 3) {
                for (aa = 0; aa < 21; aa ++)
                    for (cdn = 0; cdn < NUMS[aa]; cdn ++)
                        if (!strcmp(tri, CODONS[aa][cdn].seq)) {
                            pro[len ++] = aa;
                            bfr += log(CODONS[aa][cdn].cai);
                            goto endfor;
                        }
                
                endfor: i = 0;
            }
        }

        if (i != 0) {
            printf("Skip CDS '%.*s' due to frameshifts. \n", \
                  (int) strlen(buf) - 2, buf + 1);
            goto back;
        }

        char seq[MAX_LEN * 3];

        bfr = exp(bfr / len), aft = bfr, j = 0;

        while (j ++ < MAX_ITER) {
            for (i = 0, cai = 0; i < len; i ++) {
                ran = (float) rand() / RAND_MAX, aa = pro[i];

                for (cdn = 0; cdn < NUMS[aa]; cdn ++)
                    if (ran <= CODONS[aa][cdn].freqs) {
                        strcpy(seq + 3 * i, CODONS[aa][cdn].seq);
                        cai += log(CODONS[aa][cdn].cai);
                        break;
                    }
            }

            cai = exp(cai / len);

            coding(seq, params, strlen(seq));

            for (i = 0; i < 24; i ++)
                PyList_SetItem(vec, i, Py_BuildValue("f", params[i]));
            
            PyObject *result = PyObject_CallFunction(SVM, "O", vec);

            PyArg_Parse(result, "i", &res);

            for (i = 0; i < 4; i ++)
                if (strstr(seq, SITES[i]) != NULL)
                    break;
            
            if (i == 4 && res && cai > aft) {
                strcpy(bst, seq), aft = cai;
            }
        }

        strncpy(seq, "ATG", 3);
        fputs(buf, out);
        printf("ID: \t%.*s\nCAI:\tBefore=%.2f\tAfter=%.2f\n\n", \
              (int) strlen(buf) - 2, buf + 1, bfr, aft);

        for (i = 0, j = 0, len *= 3; i < len; i ++, j ++) {
            if (j == MAX_COL) {
                fputc('\n', out);
                j = 0;
            }

            fputc(seq[i], out);
        }

        fputc('\n', out);

        back: if (chr == '>')
            fseek(in, -1, SEEK_CUR);
    }

    Py_DECREF(pModul);
	Py_DECREF(SVM);
	Py_DECREF(vec);

    Py_Finalize();

    fclose(out);
    fclose(in);
    
    exit(0);
}