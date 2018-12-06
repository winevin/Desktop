#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>
#include <signal.h>
#include <fcntl.h>

enum
{
    ARGS = 5,
    CODE_ERR = 42,
    IN_IND = 2,
    OUT_IND = 3,
    ERR_IND = 4,
    PREM = S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP,
};

int main(int argc, char *argv[]) {
    int in, out, err;
    switch (fork()) {
    case -1:
        return 0;
    case 0:
        if (argc < ARGS) {
            return CODE_ERR;
        } else if ((in = open(argv[IN_IND], O_RDONLY)) < 0) {
            return CODE_ERR;
        } else if ((out = open(argv[OUT_IND], O_WRONLY | O_CREAT | O_APPEND, PREM)) < 0) {
            return CODE_ERR;
        } else if ((err = open(argv[ERR_IND], O_WRONLY | O_CREAT | O_TRUNC, PREM)) < 0) {
            return CODE_ERR;
        } else if (dup2(in, 0) == -1 || dup2(out, 1) == -1 || dup2(err, 2) == -1) {
            return CODE_ERR;
        } else if (close(in) || close(out) || close(err)) {
            return CODE_ERR;
        }
        execlp(argv[1], argv[1], NULL);
        return CODE_ERR;
    default:
        break;
    }
    int status;
    wait(&status);
    printf("%d\n", status);
    return 0;
}
