#include <stdio.h>
#include <dirent.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

int mycmp(const void *a, const void *b) {
    return strcasecmp(*(char **) a, *(char **) b);
}

void bypass(DIR *d, char *way, int pos) {
    int i = 0, msize = 2;
    char **dd = calloc(1, msize * sizeof(*dd));
    struct dirent *ddd;
    while ((ddd = readdir(d))) {
        if (i >= msize) {
            msize <<= 1;
            dd = realloc(dd, msize * sizeof(*dd));
        }
        if (!(ddd->d_name[0] == '.' &&
                (ddd->d_name[1] == 0 || (ddd->d_name[1] == '.' && ddd->d_name[2] == 0)))) {
            dd[i++] = strdup(ddd->d_name);
        }
    }
    closedir(d);
    qsort(dd, i, sizeof(*dd), mycmp);
    for (int j = 0; j < i; ++j) {
        char *tmp1 = way + pos, *tmp2 = dd[j];
        *tmp1++ = '/';
        while (*tmp2) {
            *tmp1++ = *tmp2++;
        }
        *tmp1 = 0;
        DIR *t = opendir(way);
        struct stat st = {0};
        if (t && lstat(way, &st) >= 0 && S_ISDIR(st.st_mode)) {
            printf("cd %s\n", dd[j]);
            bypass(t, way, tmp1 - way);
            printf("cd ..\n");
        } else if (t) {
            closedir(t);
        }
        free(dd[j]);
    }
    free(dd);
}

int main(int argc, char *argv[]) {
    char way[PATH_MAX + 1] = "";
    int k = snprintf(way, PATH_MAX, "%s", argv[1]);
    DIR *d = opendir(argv[1]);
    if (d) {
        bypass(d, way, k);
    }
    return 0;
}
