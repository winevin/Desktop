#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>

int main(void) {
    int fd[2];
    pipe(fd);
    if (!fork()) {
        if (!fork()) {
            close(fd[1]);
            long long sum = 0;
            int t, pos = 0, gg = 0;
            while (1) {
                while (pos < sizeof(t)) {
                    if ((gg = read(fd[0], &t, sizeof(t))) <= 0) {
                        break;
                    }
                    pos += gg;
                }
                if (gg <= 0) {
                    break;
                }
                pos = 0;
                sum += t;
            }
            printf("%lld\n", sum);
            close(fd[0]);
            exit(0);
        }
        close(fd[1]);
        close(fd[0]);
        while (wait(NULL) == -1);
        exit(0);
    }
    int t;
    close(fd[0]);
    int pos = 0;
    while (scanf("%d", &t) == 1) {
        while (pos < sizeof(t)) {
            pos += write(fd[1], &t, sizeof(t));
        }
        pos = 0;
    }
    close(fd[1]);
    while (wait(NULL) == -1);
    return 0;
}
