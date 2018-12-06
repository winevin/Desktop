#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "ht03-3.h"

enum
{
    MAX_RND = 62,
    CNT_NUM = 10,
    CNT_ALPHA = 26,
};

char *getline2(FILE *f) {
    int len = 2, i = 0, t;
    char *input = calloc(len + 1, sizeof(*input));
    if (input == NULL) {
        return NULL;
    }
    while ((t = getc(f)) != EOF && t != '\0') {
        input[i++] = t;
        if (i >= len) {
            len *= 2;
            char *tmp = realloc(input, len * sizeof(*input) + 1);
            if (input == NULL) {
                free(tmp);
                return NULL;
            }
            input = tmp;
        }
        if (t == '\n') {
            input[i--] = '\0';
            break;
        }
    }
    if (i == 0) {
        free(input);
        return NULL;
    }
    input[i] = '\0';
    return input;
}

int getind(RandomGenerator *rr) {
    return ((double) rr->ops->next(rr)) / MOD * MAX_RND;
}

int correct(RandomGenerator *rr, const char *pass, const char *symb) {
        const char *ptr = pass;
        while (*ptr) {
            int ind = getind(rr);
            if (*ptr != '.' && *ptr != symb[ind]) {
                break;
            }
            ++ptr;
        }
        return !*ptr;
}

int searchind(char c, const char *symb) {
    for (int i = 0; i < MAX_RND; ++i) {
        if (symb[i] == c) {
            return i;
        }
    }
    return -1;
}
int main(int argc, char *argv[]) {
    char symb[MAX_RND];
    for (int i = 0; i < CNT_NUM; ++i) {
        symb[i] = '0' + i;
    }
    for (int i = 0; i < CNT_ALPHA; ++i) {
        symb[i + CNT_NUM] = 'A' + i;
    }
    for (int i = 0; i < CNT_ALPHA; ++i) {
        symb[i + CNT_ALPHA + CNT_NUM] = 'a' + i;
    }
    char *pass = getline2(stdin);
    int seed = -1;
    int start = searchind(*pass, symb);
    RandomGenerator *rr = random_create(seed);
    for (int i = (int) ceil((double) start * MOD / MAX_RND); 
            i < (int) floor((double) (start + 1) * MOD / MAX_RND); ++i) {
        rr->now = i;
        if (correct(rr, pass + 1, symb)) {
            if (seed != -1) {
                printf("#\n");
                free(pass);
                rr->ops->destroy(rr);
                exit(0);
            }
            seed = i;
        }
    }
    if (seed == -1) {
        printf("#\n");
        free(pass);
        rr->ops->destroy(rr);
        exit(0);
    }
    char *ptr = pass + 1;
    putchar(*pass);
    rr->now = seed;
    while (*ptr) {
        putchar(symb[getind(rr)]);
        ++ptr;
    }
    putchar('\n');
    rr->ops->destroy(rr);
    free(pass);
    return 0;
}
