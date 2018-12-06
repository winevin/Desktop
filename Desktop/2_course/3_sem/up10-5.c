#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/wait.h>

int ping_pong(int max, int in, int out, int fd1, int fd2, int id) {
    if (fork()) {
        return 0;
    }
    close(fd1);
    close(fd2);
    FILE *fin = fdopen(in, "r");
    FILE *fout = fdopen(out, "w");
    setbuf(fout, NULL);
    setbuf(stdout, NULL);
    if (id == 1) {
        fprintf(fout, "1 ");
    }
    int i;
    while (fscanf(fin, "%d", &i) == 1 && i < max) {
        printf("%d %d\n", id + 1, i);
        fprintf(fout, "%d ", i + 1);
    }
    close(in);
    close(out);
    exit(0);
}

int main(int argc, char *argv[]) {
    int fd1[2], fd2[2];
    pipe(fd1);
    pipe(fd2);
    int max;
    sscanf(argv[1], "%d", &max);
    setbuf(stdout, NULL);
    ping_pong(max, fd1[0], fd2[1], fd1[1], fd2[0], 0);
    ping_pong(max, fd2[0], fd1[1], fd1[0], fd2[1], 1);
    close(fd1[0]);
    close(fd1[1]);
    close(fd2[0]);
    close(fd2[1]);
    wait(NULL);
    wait(NULL);
    printf("Done\n");
    return 0;
}
