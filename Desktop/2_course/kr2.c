#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/wait.h>

int myexec(char *file, int fd[], int id) {
    if (fork()) {
        return 0;
    }
    close(fd[0]);
    unsigned ans = 0, f = open(file, O_RDONLY);
    dup2(f, 0);
    dup2(fd[1], 1);
    close(f);
    while (scanf("%d", &f) == 1) {
        ans += f;
    }
    printf("%d %d ", id, ans);
    fflush(stdout);
    exit(0);
}

int main(int argc, char *argv[]) {
    int fd[2];
    pipe(fd);
    myexec(argv[1], fd, 0);
    myexec(argv[2], fd, 1);
    dup2(fd[0], 0);
    close(fd[0]);
    close(fd[1]);
    wait(NULL);
    wait(NULL);
    unsigned id[2], ans[2];
    scanf("%d%d", id, ans);
    scanf("%d%d", id + 1, ans + 1);
    printf("%d\n%d\n%d\n", ans[1 - id[0]], ans[1 - id[1]], ans[0] + ans[1]);
    return 0;
}
