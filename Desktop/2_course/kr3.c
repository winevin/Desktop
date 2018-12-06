#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>
#include <string.h>

int myexec(char **argv) {
    pid_t pid = fork();
    if (pid == -1) {
        exit(1);
    } else if (!pid) {
        execvp(argv[0], argv);
        exit(1);
    }
    int status;
    wait(&status);
    return WIFEXITED(status) && !WEXITSTATUS(status);
}

int main(int argc, char *argv[]) {
    char *sep = argv[1];
    int *start = calloc(argc, sizeof(*start)); 
    int size = 0, flag = 1, flagA = 1;
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(sep, argv[i])) {
            if (flag) {
                start[size++] = i + 1;
                flag = 0;
            }
            argv[i] = NULL;
        } else {
            flagA = 0;
            flag = 1;
        }
    }
    if (flagA) {
        free(start);
        return 0;
    }
    for (int j = 0; j < size; ++j) {
        if (!myexec(argv + start[j])) {
            free(start);
            return 1;
        }
    }
    free(start);
    return 0;
}
