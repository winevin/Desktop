#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/wait.h>

int max, id, *pid;
FILE *fin, *fout;

void err(int e) {
    if (pid == NULL) {
        perror("pid is NULL");
        exit(1);
    } else if (fin == NULL) {
        perror("fin is NULL");
        exit(1);
    } else if (fout == NULL) {
        perror("fout is NULL");
        exit(1);
    } else {
        perror("Just exit");
        exit(1);
    }
}

void mysig(int sig) {
    signal(sig, mysig);
    int i;
    if (fscanf(fin, "%d", &i) != 1 || i >= max) {
    err(1);
        fclose(fin);
        fclose(fout);
        exit(0);
    }
    printf("%d %d\n", id + 1, i);
    fprintf(fout, "%d ", i + 1);
    kill(pid[id], SIGUSR1);
}

int ping_pong(int fd[]) {
    int tmppid = fork();
    if (tmppid) {
        return tmppid;
    }
    fin = fdopen(fd[0], "r");
    fout = fdopen(fd[1], "w");
    setbuf(fout, NULL);
    setbuf(stdout, NULL);
    if (id == 1) {
        fprintf(fout, "1 ");
    }
    while (1) {
        pause();
    }
}

int main(int argc, char *argv[]) {
    int fd[2];
    pid = calloc(2, sizeof(*pid));
    signal(SIGUSR1, mysig);
    pipe(fd);
    sscanf(argv[1], "%d", &max);
    setbuf(stdout, NULL);
    pid[1] = ping_pong(fd);
    id = 1;
    pid[0] = ping_pong(fd);

    kill(pid[1], SIGUSR1);

    close(fd[0]);
    close(fd[1]);
    wait(NULL);
    wait(NULL);
    free(pid);
    printf("Done\n");
    return 0;
}
