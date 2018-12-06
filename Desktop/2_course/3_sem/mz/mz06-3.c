#include <stdlib.h>
#include <string.h>

int parse_rwx_permissions(const char *str) {
    const char temp[] = "rwxrwxrwx";
    int ans = 0;
    if (str == NULL || strlen(str) != sizeof(temp) - 1) {
        return -1;
    }
    for (int i = 0; i < sizeof(temp) - 1; ++i) {
        ans <<= 1;
        if (str[i] == temp[i]) {
            ++ans;
        } else if (str[i] != '-') {
            return -1;
        }
    }
    return ans;
}
