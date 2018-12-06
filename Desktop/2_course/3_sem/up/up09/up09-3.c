#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>

int proc(char *cmd) {
    if (!fork()) {
        execlp(cmd, cmd, NULL);
        exit(1);
    }
    int status;
    wait(&status);
    return WIFEXITED(status) && !WEXITSTATUS(status);
}

int main(int argc, char *argv[]) {
    return !((proc(argv[1]) || proc(argv[2])) && proc(argv[3]));
}
