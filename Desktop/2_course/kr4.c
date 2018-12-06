#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/wait.h>
#include <math.h>

enum
{
    IDX_D = 1,
    IDX_P,
    IDX_O,
    IDX_X,
    IDX_N,
    ACCESS = S_IWUSR | S_IRUSR | S_IWGRP | S_IRGRP | S_IWOTH | S_IROTH,
    SETTINGS = O_WRONLY,
    SHIFT = 18,
};

int myfork(int d, int p, char *o, int x, int n, int id) {
    if (fork()) {
        return 0;
    }
    int fd = open(o, SETTINGS);
    FILE *f = fdopen(fd, "w+");
    int a = x + id;
    for (int i = id; i < n; i += p) {
        lseek(fd, i * SHIFT, SEEK_SET);
        fprintf(f, "%+10e\n", cos(a));
        a = x + i * d;
    }
    fclose(f);
    return 0;
}

int main(int argc, char *argv[]) {
    int d, p, x, n;
    sscanf(argv[IDX_D], "%d", &d);
    sscanf(argv[IDX_P], "%d", &p);
    sscanf(argv[IDX_X], "%d", &x);
    sscanf(argv[IDX_N], "%d", &n);
    close(open(argv[IDX_O], SETTINGS | O_CREAT | O_TRUNC, ACCESS));
    for (int i = 0; i < p; ++i) {
        myfork(d, p, argv[IDX_O], x, n, i);
    }
    for (int i = 0; i < p; ++i) {
        wait(NULL);
    }
    return 0;
}
