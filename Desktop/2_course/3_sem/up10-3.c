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
    FILE_B,
    ACCESS = S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP,
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
        int file1 = open(argv[FILE_A], O_RDONLY);
        close(fd[0]);
        if (myexec(argv[CMD_A], file1, fd[1])) {
            myexec(argv[CMD_B], 0, fd[1]);
        }
        close(file1);
        close(fd[1]);
        exit(0);
    }
    if (!fork()) {
        int file2 = open(argv[FILE_B], FLAG_WRITE, ACCESS);
        close(fd[1]);
        myexec(argv[CMD_C], fd[0], file2);
        close(file2);
        close(fd[0]);
        exit(0);
    }
    close(fd[0]);
    close(fd[1]);
    wait(NULL);
    wait(NULL);
    return 0;
}
