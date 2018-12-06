#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
    int cnt, low, high, seed;
    sscanf(argv[1], "%d", &cnt);
    sscanf(argv[2], "%d", &low);
    sscanf(argv[3], "%d", &high);
    sscanf(argv[4], "%d", &seed);
    srand(seed);
    for (int i = 0; i < cnt; ++i) {
        printf("%d\n", (int) ((double) rand() / ((double) RAND_MAX + 1) * (high - low)) + low);
    }
    return 0;
}
