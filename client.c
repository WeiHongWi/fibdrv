#include <fcntl.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#define FIB_DEV "/dev/fibonacci"

static inline long long u_ntime()
{
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return ts.tv_sec * 1e9 + ts.tv_nsec;
}

int main()
{
    char buf[100000];
    char write_buf[] = "testing writing";
    int offset = 2000; /* TODO: try test something bigger than the limit */


    FILE *data = fopen("datav2.txt", "w");
    int fd = open(FIB_DEV, O_RDWR);

    if (fd < 0) {
        perror("Failed to open character device");
        exit(1);
    }
    write(fd, write_buf, strlen(write_buf));
    for (int i = 0; i <= offset; i++) {
        lseek(fd, i, SEEK_SET);
        long long start = u_ntime();
        long long sz = read(fd, buf, 1);
        long long utime = u_ntime() - start;
        long long ktime = write(fd, write_buf, strlen(write_buf));

        // fprintf(data, "%d %lld %lld %lld\n", i,utime,ktime,utime - ktime);
        fprintf(data, "%d %lld \n", i, ktime);
        printf("Reading from " FIB_DEV
               " at offset %d, returned the sequence "
               "%s.\n",
               i, buf);
    }

    close(fd);
    fclose(data);
    return 0;
}
