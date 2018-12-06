#include <stdio.h>
#include <stdlib.h>

enum { PERCENT = 100 };

int main(void) {
    int n, m, seed;
    scanf("%d", &n);
    int *qsi = calloc(n, sizeof(*qsi));
    int *ver = calloc(n, sizeof(*ver));
    if (!qsi || !ver) {
        return 1;
    }
    for (int i = 0; i < n; ++i) {
        scanf("%d%d", qsi + i, ver + i);
        if (ver[i] == 0) {
            --i;
            --n;
            continue;
        }
        if (i) {
            ver[i] += ver[i - 1];
        }
    }
    scanf("%d%d", &m, &seed);
    srand(seed);
    for (int i = 0; i < m; ++i) {
        double val = ((double) rand() / (((double) RAND_MAX) + 1) * PERCENT);
        for (int j = 0; j < n; ++j) {
            if ((double)ver[j] >= val) {
                printf("%d\n", qsi[j]);
                break;
            }
        }
    }
    free(qsi);
    free(ver);
    return 0;
}
