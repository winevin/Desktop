#include <stdio.h>
#include <stdlib.h>
#include <sys/wait.h>
#include <unistd.h>
#include <fcntl.h>

enum
{
    CMD_A = 1,
    CMD_B,
    CMD_C,
    CMD_D,
    FILE_A,
    ACCESS = S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP,
    FLAG_WRITE = O_WRONLY | O_TRUNC | O_CREAT,
};

int myexec(char *cmd, int in, int out) {
    if (!fork()) {
        if (in != 0) {
            dup2(in, 0);
            close(in);
        }
        dup2(out, 1);
        close(out);
        execlp(cmd, cmd, NULL);
        exit(EXIT_FAILURE);
    }
    int status;
    wait(&status);
    return WIFEXITED(status) && !WEXITSTATUS(status);
}

int main(int argc, char *argv[]) {
    int fd[2];
    pipe(fd);
    if (!fork()) {
        close(fd[0]);
        if (myexec(argv[CMD_A], 0, fd[1])) {
            myexec(argv[CMD_B], 0, fd[1]);
        }
        close(fd[1]);
        exit(0);
    }
    if (!fork()) {
        close(fd[1]);
        int f = open(argv[FILE_A], FLAG_WRITE, ACCESS);
        myexec(argv[CMD_C], fd[0], f);
        myexec(argv[CMD_D], fd[0], f);
        close(f);
        close(fd[0]);
    }
    close(fd[0]);
    close(fd[1]);
    wait(NULL);
    wait(NULL);
    return 0;
}
