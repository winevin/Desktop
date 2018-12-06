#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>

enum
{
    CNT_NUM = 10,
    LEN_MAX = PATH_MAX + 2 //for \n && \r
};

int main(void) {
    char in[LEN_MAX];
    if (fgets(in, sizeof(in), stdin) == NULL) {
        in[0] = '\0';
    }
    int len = strlen(in);
    while (len > 0 && (in[len - 1] == '\n' || in[len - 1] == '\r')) {
        in[--len] = '\0';
    }
    unsigned long long cnt[CNT_NUM] = {0};
    if (len > 0) {
        FILE *f = fopen(in, "r");
        if (f != NULL) {
            int t;
            while ((t = getc_unlocked(f)) != EOF) {
                if (t >= '0' && t <= '9') {
                    cnt[t - '0']++;
                }
            }
            fclose(f);
        }
    }
    for (int j = 0; j < CNT_NUM; j++) {
        printf("%d %llu\n", j, cnt[j]);
    }
    return 0;
}
