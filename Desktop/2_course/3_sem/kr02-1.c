#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>

static inline int myexec(char *cmd) {
    if (!fork()) {
        execlp(cmd, cmd, NULL);
        exit(1);
    }
    return 0;
}

static inline mywait(void) {
    int status;
    while (wait(&status) == -1);
    return WIFEXITED(status) && !WEXITSTATUS(status);
}

int main(int argc, char *argv[]) {
    for (int i = 1; i < argc; ++i) {
        myexec(argv[i]);
        if (!mywait()) {
            return 1;
        }
    }
    return 0;
}
