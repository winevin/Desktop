void normalize_path(char *buf) {
    char *p1 = buf, *p2 = buf;
    while (*++p1) {
        if (*p1 != '/' || *p2 != '/') {
            *++p2 = *p1;
        }
    }
    *++p2 = '\0';
}
