#include <stdio.h>
#include <stdlib.h>

char * getline2(FILE *f) {
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
