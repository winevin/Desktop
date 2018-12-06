#include <stdio.h>
#include <stdlib.h>

enum { MAX_LEN = 16384 };

typedef struct Node
{
    double first;
    int second;
    int index;
} Pair;

int cmp(const void *a, const void *b) {
    if (((Pair *)a)->first > ((Pair *)b)->first) {
        return 1;
    }
    if (((Pair *)a)->first < ((Pair *)b)->first) {
        return -1;
    }
    return ((Pair *)a)->index - ((Pair *)b)->index;
}

int main(void) {
    int n = 0;
    Pair *in = calloc(MAX_LEN, sizeof(*in));
    if (in == NULL) {
        return 2;
    }
    while (scanf("%lf%d", &(in[n].first), &(in[n].second)) == 2 && n < MAX_LEN) {
        in[n].index = n;
        ++n;
    }
    qsort(in, n, sizeof(*in), cmp);
    for (int i = 0; i < n; i++) {
        printf("%.10g %d\n", in[i].first, in[i].second);
    }
    free(in);
    return 0;
}
