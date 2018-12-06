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
    FILE_A,
    ACCESS = S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH,
    FLAG_WRITE = O_WRONLY | O_APPEND | O_CREAT,
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
        myexec(argv[CMD_A], 0, fd[1]);
        close(fd[1]);
        exit(0);
    }
    if (!fork()) {
        int file = open(argv[FILE_A], FLAG_WRITE, ACCESS);
        close(fd[1]);
        if (myexec(argv[CMD_B], fd[0], file)) {
            myexec(argv[CMD_C], fd[0], file);
        }
        close(file);
        close(fd[0]);
        exit(0);
    }
    close(fd[0]);
    close(fd[1]);
    wait(NULL);
    wait(NULL);
    return 0;
}
