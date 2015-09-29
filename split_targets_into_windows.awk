# splits targets >= 1000 bp long
# into equally sized windows
# of a size within [500, 1000) bp

{
    len = $3 - $2;
    if (len < 1000) {
        printf "%s\t%d\t%d\t%s:%d-%d\t%s:%d-%d\n",
                $1, $2, $3, $1, $2, $3, $1, $2, $3;
        next;
    }

    n = int(len / 500);
    spacing = len / n;

    for (i = 0; i < n; i++) {
        start = int($2 + spacing * i);
        end   = int($2 + spacing * (i+1));
        printf "%s\t%d\t%d\t%s:%d-%d\t%s:%d-%d\n",
                $1, start, end, $1, start, end, $1, $2, $3;
    }
}

