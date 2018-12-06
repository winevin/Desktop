#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

enum
{
    MAX_SEP = 0x20,
    FOUR_BIT = 4,
};

static inline unsigned char *getline2(FILE *f, int start_size) {
    int len = 2 + start_size, i = 0, t = 0;
    unsigned char *input = calloc(len + 1, sizeof(*input));
    if (input == NULL) {
        return NULL;
    }
    while (t <= MAX_SEP) {
        t = getc_unlocked(f);
        if (t == -1) {
            free(input);
            return NULL;
        }
    }
    if (t > MAX_SEP) {
        input[i++] = t;
    }
    while ((t = getc_unlocked(f)) != EOF && t > MAX_SEP) {
        input[i++] = t;
        if (i >= len) {
            len <<= 1;
            unsigned char *tmp = realloc(input, len * sizeof(*tmp) + 1);
            if (tmp == NULL) {
                free(input);
                return NULL;
            }
            input = tmp;
        }
    }
    if (i == 0) {
        free(input);
        return NULL;
    }
    input[i] = '\0';
    return input;
}

static inline int utf_byte(const unsigned char c) {
    int c1 = c >> (CHAR_BIT - FOUR_BIT);
    int ans = 0;
    for (int i = 0; i < FOUR_BIT; ++i) {
        if (c1 & 1) {
            ++ans;
        } else {
            ans = 0;
        }
        c1 >>= 1;
    }
    return ans ? ans : ans + 1;
}

static inline int strlen2(const unsigned char *s) {
    int len = 0;
    while (*s) {
        s += utf_byte(*s);
        ++len;
    }
    return len;
}

int main(void) {
    unsigned char *s, *ans = NULL;
    int max_len = 0, len;
    while ((s = getline2(stdin, max_len)) != NULL) {
        len = strlen2(s);
        if (len > max_len) {
            max_len = len;
            free(ans);
            ans = s;
        } else {
            free(s);
        }
    }
    printf("%d\n", max_len);
    if (max_len > 0) {
        printf("%s", ans);
    }
    putchar('\n');
    free(ans);
    return 0;
}
