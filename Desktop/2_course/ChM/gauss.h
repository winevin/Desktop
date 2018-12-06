#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "svd.h"

#define IDX(x, y) (x * my + y)
#define DEBUG 0X0

#define EPS 1e-10

int print(double *matrix, int mx, int my);
int DiagNorm(double *matrix, int mx, int my, int pos); //0x1
int SubDiag(double *matrix, int mx, int my, int pos, int start, int end); //0x2
double *CopyMatrix(double *matrix, int n);
int Swap2Lines(double *matrix, int mx, int my, int a, int b);
int Swap2Column(double *matrix, int mx, int my, int a, int b);
double * Gauss(double *matrix, double *ans, int mx, int my); //0x8
double * GaussWithMain(double *matrix, double *ans, int mx, int my); //0x10

int Print(double *matrix, int mx, int my) {
    printf("\n------------------------------------\n\n");
    for (int i = 0; i < mx; ++i) {
        for (int j = 0; j < my; ++j) {
            printf("%10.06lf ", matrix[IDX(i, j)]);
        }
        printf("\n");
    }
    printf("\n++++++++++++++++++++++++++++++++++++\n");
    return 0;
}

int DiagNorm(double *matrix, int mx, int my, int pos) {
    if (pos >= mx || pos >= my) {
        printf("DiagNorm(%p, %d, %d, %d)", matrix, mx, my, pos);
        perror("in DiagNorm pos >= mx or pos >= my");
        abort();
    }
    double diag = matrix[IDX(pos, pos)];
    if (fabs(diag) < 0.00001) {
        perror("in DiagNorm diag elem == 0");
        abort();
    }
    for (int i = 0; i < my; ++i) {
        matrix[IDX(pos, i)] /= diag;
    }
    return DEBUG & 0x1 && printf("\nvvv DiagNorm vvv") && Print(matrix, mx, my);
}

int SubDiag(double *matrix, int mx, int my, int pos, int start, int end) {
    for (int i = start; i < end; ++i) {
        double cur = matrix[IDX(i, pos)];
        for (int j = IDX(i, 0), jpos = IDX(pos, 0); j < IDX(i, my); ++j, ++jpos) {
            matrix[j] -= matrix[jpos] * cur;
        }
    }
    return DEBUG & 0x2 && printf("\nvvv SubDiag vvv") && Print(matrix, mx, my);
}

double *CopyMatrix(double *matrix, int mx) {
    int my = mx * 2 + 1;
    double *m = calloc(my * mx, sizeof(*m));
    for (int i = 0; i < mx; ++i) {
        for (int j = 0; j < mx + 1; ++j) {
            m[IDX(i, j)] = matrix[i * (mx + 1) + j];
        }
        m[IDX(i, i + mx + 1)] = 1;
    }
    return m;
}

int Swap2Lines(double *matrix, int mx, int my, int a, int b) {
    for (int i = IDX(a, 0), j = IDX(b, 0); i < IDX(a, my); ++i, ++j) {
        double tmp = matrix[i];
        matrix[i] = matrix[j];
        matrix[j] = tmp;
    }
    return 0;
}

int Swap2Column(double *matrix, int mx, int my, int a, int b) {
    for (int i = IDX(0, a), j = IDX(0, b); i < IDX(mx, a); i += my, j += my) {
        double tmp = matrix[i];
        matrix[i] = matrix[j];
        matrix[j] = tmp;
    }
    return 0;
}

int GetAns(double *matrix, double *ans, int mx, int my) {
    for (int i = mx - 1; i >= 0; --i) {
        ans[i] = matrix[IDX(i, mx)];
        for (int j = i + 1; j < mx; ++j) {
            ans[i] -= matrix[IDX(i, j)] * ans[j];
        }
    }
    return 0;
}

double * Invert(double *matrix, int mx, int my) {
    double *m = calloc(mx * mx, sizeof(*m));
    for (int i = 0; i < mx; ++i) {
        for (int j = 0; j < mx; ++j) {
            m[i * mx + j] = matrix[IDX(i, j + mx + 1)];
        }
    }
    free(matrix);
    return m;
}

double * Gauss(double *matrix, double *ans, int mx, int my) {
    matrix = CopyMatrix(matrix, mx);
    for (int i = 0; i < mx; ++i) {
        int j;
        for (j = i; j < mx; ++j) {
            if (fabs(matrix[IDX(j, i)]) > EPS) {
                break;
            }
        }
        if (j != i) {
            Swap2Lines(matrix, mx, my, i, j);
        }
        DiagNorm(matrix, mx, my, i);
        SubDiag(matrix, mx, my, i, i + 1, mx);
    }
    GetAns(matrix, ans, mx, my);
    for (int i = mx - 1; i >= 0; --i) {
        SubDiag(matrix, mx, my, i, 0, i);
    }
    return !(DEBUG & 0x8) || printf("\nvvv Gauss vvv") == 0 || Print(matrix, mx, my), Invert(matrix, mx, my);
}

double * GaussWithMain(double *matrix, double *ans, int mx, int my) {
    matrix = CopyMatrix(matrix, mx);
    double answer[mx];
    int map[mx];
    for (int i = 0; i < mx; ++i) {
        map[i] = i;
    }
    for (int i = 0; i < mx; ++i) {
        int mj = i;
        for (int j = i; j < mx; ++j) {
            if (fabs(matrix[IDX(i, j)]) > fabs(matrix[IDX(i, mj)])) {
                mj = j;
            }
        }
        if (mj != i) {
            map[i] = mj;
            map[mj] = i;
            Swap2Column(matrix, mx, my, i, mj);
            if (mx * 2 + 1 == my) {
                Swap2Column(matrix, mx, my, my - mx +  i, my - mx + mj);
            }
        }
        DiagNorm(matrix, mx, my, i);
        SubDiag(matrix, mx, my, i, i + 1, mx);
    }
    GetAns(matrix, answer, mx, my);
    for (int i = 0; i < mx; ++i) {
        ans[i] = answer[map[i]];
    }
    for (int i = mx - 1; i >= 0; --i) {
        SubDiag(matrix, mx, my, i, 0, i);
    }
    return !(DEBUG & 0x10) || printf("\nvvv GaussWithMain vvv") == 0 || Print(matrix, mx, my), Invert(matrix, mx, my);
}

double NumOfObus(double *matrix, int mx) {
    float w[mx], v[mx][mx], m[mx][mx];
    double dw[mx];
    for (int i = 0; i < mx; ++i)
        for (int j = 0; j < mx; ++j)
            m[i][j] = matrix[mx * i + j];
    dsvd(m, mx, mx, w, v);
    for (int i = 0; i < mx; ++i)
        dw[i] = w[i];
    Print(dw, mx, 1);
    return 0;
}
