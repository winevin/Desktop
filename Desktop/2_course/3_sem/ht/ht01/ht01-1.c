int bitcount(STYPE value) {
    UTYPE num = value;
    int cnt = 0;
    while (num) {
        cnt += num & 1;
        num >>= 1;
    }
    return cnt;
}
