#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>

enum { TERMSIG_START = 128 };

int mysys(const char *str) {
    pid_t pid = fork();
    if (pid == -1) {
        return -1;
    } else if (!pid) {
        execlp("/bin/sh", "sh", "-c", str, NULL);
        _exit(TERMSIG_START - 1);
    }
    int status;
    waitpid(pid, &status, 0);
    if (WIFEXITED(status)) {
        return WEXITSTATUS(status);
    } else if (WIFSIGNALED(status)) {
        return TERMSIG_START + WTERMSIG(status);
    }
    abort();
}
