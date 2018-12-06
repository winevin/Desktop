#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
    DIR *d = opendir(argv[1]);
    if (d == NULL) {
        return 1;
    }
    struct dirent *dd;
    uid_t uid = getuid();
    struct stat st;
    long long ans = 0;
    while ((dd = readdir(d))) {
        char name[PATH_MAX + 1];
        if (snprintf(name, PATH_MAX + 1, "%s/%s", argv[1], dd->d_name) < PATH_MAX
                && stat(name, &st) == 0
                && dd->d_name[0] >= 'A' && dd->d_name[0] <= 'Z'
                && S_ISREG(st.st_mode) && st.st_uid == uid) {
            ans += st.st_size;
        }
    }
    printf("%lld\n", ans);
    closedir(d);
    return 0;
}
