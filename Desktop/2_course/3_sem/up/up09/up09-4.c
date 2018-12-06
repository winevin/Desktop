#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <limits.h>
#include <time.h>
#include <unistd.h>

enum { SIZE = 33 };

int main(int argc, char *argv[]) {
    char tmp[] = "/tmp";
    char *tmp_path = getenv("XDG_RUNTIME_DIR");
    if (tmp_path == NULL) {
        tmp_path = getenv("TMPDIR");
        if (tmp_path == NULL) {
            tmp_path = tmp;
        }
    }
    char path[PATH_MAX + 1] = {0}, name[SIZE + 1] = {0};
    srand(time(NULL));
    snprintf(name, SIZE, "%016lx%016x", time(NULL), rand());
    snprintf(path, PATH_MAX, "%s/%s", tmp_path, name);
    int fd = open(path, O_WRONLY | O_TRUNC | O_CREAT, S_IRWXU);
    FILE *f = fdopen(fd, "w");
    fprintf(f, "#!/usr/bin/env python3\nprint(");
    for (int i = 1; i < argc - 1; ++i) {
        fprintf(f, "%s*", argv[i]);
    }
    fprintf(f, "%s)\n", argv[argc - 1]);
    fclose(f);
    execlp(path, name, NULL);
    return 1;
}
