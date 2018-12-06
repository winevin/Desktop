#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/wait.h>

enum
{
    IDX_N = 1,
    IDX_F,
    IDX_A,
    IDX_D,
    IDX_K,
    ACCESS = S_IWUSR | S_IRUSR | S_IWGRP | S_IRGRP | S_IWOTH | S_IROTH,
    SETTINGS = O_WRONLY,
};

int myfork(int n, char *f, int a, int d, int k, int id) {
    if (fork()) {
        return 0;
    }
    int fd = open(f, SETTINGS);
    lseek(fd, id * sizeof(a), SEEK_SET);
    a += d * id;
    for (int i = 0; i < k; ++i) {
        write(fd, &a, sizeof(a));
        a += d * n;
        lseek(fd, (n - 1) * sizeof(a), SEEK_CUR);
    }
    close(fd);
    return 0;
}

int main(int argc, char *argv[]) {
    int n, a, d, k;
    sscanf(argv[IDX_N], "%d", &n);
    sscanf(argv[IDX_A], "%d", &a);
    sscanf(argv[IDX_D], "%d", &d);
    sscanf(argv[IDX_K], "%d", &k);
    close(open(argv[IDX_F], SETTINGS | O_CREAT | O_TRUNC, ACCESS));
    for (int i = 0; i < n; ++i) {
        myfork(n, argv[IDX_F], a, d, k, i);
    }
    for (int i = 0; i < n; ++i) {
        wait(NULL);
    }
    return 0;
}
