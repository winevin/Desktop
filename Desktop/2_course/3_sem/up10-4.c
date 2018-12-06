#include <unistd.h>
#include <stdlib.h>
#include <sys/wait.h>
#include <signal.h>

int killall(pid_t *pid, int i) {
    for (int j = 1; j < i; ++j) {
        if (pid[i]) {
            kill(pid[i], SIGKILL);
        }
    }
    free(pid);
    exit(1);    
}

inline pid_t myexec(char *cmd, int in, int out, int fd) {
    int pid;
    if (!(pid = fork())) {
        if (in != 0) {
            dup2(in, 0);
            close(in);
        }
        if (out != 1) {
            dup2(out, 1);
            close(out);
        }
        if (fd != 0) {
            close(fd);
        }
        execlp(cmd, cmd, NULL);
        exit(1);
    }
    return pid;
}

int main(int argc, char *argv[]) {
    int prev_fd0;
    int next_fd[2] = {0, 1};
    pid_t *pid = calloc(argc + 1, sizeof(*pid));
    if (pid == NULL) {
        exit(1);
    }
    for (int i = 1; i < argc; ++i) {
        prev_fd0 = next_fd[0];
        if (i == argc - 1) {
            next_fd[0] = 0;
            next_fd[1] = 1;
        } else {
            if (pipe(next_fd)) {
                killall(pid, i);
            }
        }

        if ((pid[i] = myexec(argv[i], prev_fd0, next_fd[1], next_fd[0])) == -1) {
            killall(pid, i);
        }
        if ((next_fd[1] != 1 && close(next_fd[1])) || (prev_fd0 != 0 && close(prev_fd0))) {
            killall(pid, i + 1);
        }
    }
    for (int i = 1; i < argc; ++i) {
        if (wait(NULL) == -1) {
            killall(pid, argc);
        }
    }
    free(pid);
    return 0;
}
