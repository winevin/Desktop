#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>

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

int myexec(char *file) {
    FILE *f = fopen(file, "r");
    char *cmd = getline2(f);
    fclose(f);
    if (!fork()) {
        execlp(cmd, cmd, NULL);
        exit(1);
    }
    free(cmd);
    return 0;
}

int mywait(void) {
    int status;
    while (wait(&status) == -1);
    return WIFEXITED(status) && !WEXITSTATUS(status);
}

int main(int argc, char *argv[]) {
    long long n, cnt = 0;
    sscanf(argv[1], "%lld", &n);
    n = n ? n : 1;
    n = n < argc - 2 ? n : argc - 2;
    for (int i = 1; i < n + 1; ++i) {
        myexec(argv[i + 1]);
    }
    for (int i = 1; i < n + 1; ++i) {
        cnt += mywait();
    }
    for (int i = n + 1; i < argc - 1; ++i) {
        myexec(argv[i + 1]);
        cnt += mywait();
    }
    printf("%lld\n", cnt);
    return 0;
}
