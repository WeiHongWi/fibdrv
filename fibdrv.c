#include <linux/cdev.h>
#include <linux/device.h>
#include <linux/fs.h>
#include <linux/init.h>
#include <linux/kdev_t.h>
#include <linux/kernel.h>
#include <linux/module.h>
#include <linux/mutex.h>
#include <linux/slab.h>
#include "bn_kernel.c"


MODULE_LICENSE("Dual MIT/GPL");
MODULE_AUTHOR("National Cheng Kung University, Taiwan");
MODULE_DESCRIPTION("Fibonacci engine driver");
MODULE_VERSION("0.1");

#define DEV_FIBONACCI_NAME "fibonacci"



#define MAX_LENGTH 10000

static dev_t fib_dev = 0;
static struct cdev *fib_cdev;
static struct class *fib_class;
static DEFINE_MUTEX(fib_mutex);
static ktime_t kt;


void reverse_str(char *a)
{
    size_t len = strlen(a);
    for (int i = 0; i < len / 2; ++i) {
        char temp = a[i];
        a[i] = a[len - i - 1];
        a[len - i - 1] = temp;
    }
}

typedef struct Bignum {
    char element[128];
} bn_str;


static void string_number_add(char *a, char *b, char *out)
{
    reverse_str(a);
    reverse_str(b);
    int i;
    int carry = 0, sum;
    size_t size_a = strlen(a);
    size_t size_b = strlen(b);
    for (i = 0; i < size_b; i++) {
        sum = (a[i] - '0') + (b[i] - '0') + carry;
        out[i] = '0' + sum % 10;
        carry = sum / 10;
    }
    for (i = size_b; i < size_a; i++) {
        sum = (a[i] - '0') + carry;
        out[i] = '0' + sum % 10;
        carry = sum / 10;
    }
    if (carry)
        out[i++] = '0' + carry;

    out[i] = '\0';
    reverse_str(a);
    reverse_str(b);
    reverse_str(out);
}
/*
static long long pow1(long long base, int expo)
{
    long long num = 1;
    for (int i = 0; i < expo; ++i) {
        num = num * base;
    }
    return num;
}
static long long fib_sequence(long long k)
{
    long long a = 0, b = 1;
    int lz_bit = __builtin_clzll(k);
    long long current_digit = 1 << (64 - lz_bit - 1);

    for (int i = 64 - lz_bit; i >= 1; i--) {
        long long t1 = a * (2 * b - a);
        long long t2 = pow1(a, 2) + pow1(b, 2);
        a = t1;
        b = t2;
        if (k & current_digit) {
            t1 = a + b;
            a = b;
            b = t1;
        }
        current_digit >>= 1;
    }
    return a;
}
*/

void fib_sequence_bn(long long k, bn *dest)
{
    bn_resize(dest, 1);
    if (k < 2) {  // Fib(0) = 0, Fib(1) = 1
        dest->number[0] = k;
        return;
    }

    bn *a = bn_alloc(1);
    bn *b = bn_alloc(1);
    dest->number[0] = 1;

    for (unsigned int i = 1; i < k; i++) {
        bn_swap(b, dest);
        bn_add(a, b, dest);
        bn_swap(a, b);
    }
    bn_free(a);
    bn_free(b);
}
void fib_sequence_fastd(long long k, bn *dest)
{
    bn_resize(dest, 1);
    if (k < 2) {  // Fib(0) = 0, Fib(1) = 1
        dest->number[0] = k;
        return;
    }

    bn *f1 = bn_alloc(1);
    bn *f2 = dest;
    f1->number[0] = 0;
    f2->number[0] = 1;
    bn *k1 = bn_alloc(1);
    bn *k2 = bn_alloc(1);

    for (unsigned int i = 1U << (30 - __builtin_clz(k)); i; i >>= 1) {
        /* F(2k) = F(k) * [ 2 * F(k-1) + F(k) ] */
        bn_add(f1, f1, k1);
        bn_add(k1, f2, k1);
        bn_mult(k1, f2, k2);
        // F(2k-1) = F(k)^2 + F(k-1)^2
        bn_mult(f2, f2, k1);
        bn_swap(f2, k2);
        bn_mult(f1, f1, k2);
        bn_add(k1, k2, f1);
        // bn_swap(f1,k1);
        if (k & i) {
            bn_swap(f1, f2);
            bn_add(f1, f2, f2);
        }
    }
    bn_free(f1);
    bn_free(k1);
    bn_free(k2);
}
static long long fib_sequence_str(long long k, char *buf)
{
    bn_str *f = kmalloc(sizeof(bn_str) * (k + 2), GFP_KERNEL);

    strncpy(f[0].element, "0", 1);
    f[0].element[1] = '\0';
    strncpy(f[1].element, "1", 1);
    f[1].element[1] = '\0';
    for (int i = 2; i <= k; ++i) {
        string_number_add(f[i - 1].element, f[i - 2].element, f[i].element);
    }
    unsigned long ret_size = strlen(f[k].element) + 1;
    unsigned long long test = __copy_to_user(buf, f[k].element, ret_size);
    if (test) {
        printk("The copy from kernel to user is fail.");
        return -1;
    }
    return ret_size;
}
// cppcheck-suppress unusedFunction
static long long fib_time_proxy(long long k, char *buf, int ctrl)
{
    long long result = 0;
    kt = ktime_get();
    switch (ctrl) {
    case 1:
        // result = fib_sequence_str(k, buf);
        break;
    case 2:
        // fib_sequence_bn(k, buf);
        break;
    case 3:
        // result = fib_sequence_fastd(k, buf);
        break;
    }
    kt = ktime_sub(ktime_get(), kt);

    return result;
}
static int fib_open(struct inode *inode, struct file *file)
{
    if (!mutex_trylock(&fib_mutex)) {
        printk(KERN_ALERT "fibdrv is in use");
        return -EBUSY;
    }
    return 0;
}

static int fib_release(struct inode *inode, struct file *file)
{
    mutex_unlock(&fib_mutex);
    return 0;
}

/* calculate the fibonacci number at given offset */
static ssize_t fib_read(struct file *file,
                        char *buf,
                        size_t size,
                        loff_t *offset)
{
    bool doubling = true;
    bn *res = bn_alloc(1);
    if (doubling) {
        kt = ktime_get();
        fib_sequence_fastd(*offset, res);
        kt = ktime_sub(ktime_get(), kt);
    } else {
        kt = ktime_get();
        fib_sequence_bn(*offset, res);
        kt = ktime_sub(ktime_get(), kt);
    }
    char *p = bn_to_string(res);
    size_t sz = strlen(p) + 1;
    access_ok(buf, size);
    __copy_to_user(buf, p, sz);
    bn_free(res);

    return 1;

    // return (ssize_t) fib_time_proxy(*offset, buf, ctrl);
    //  return fib_sequence(*offset,buf);
}
/* write operation is skipped */
static ssize_t fib_write(struct file *file,
                         const char *buf,
                         size_t size,
                         loff_t *offset)
{
    return ktime_to_ns(kt);
}

static loff_t fib_device_lseek(struct file *file, loff_t offset, int orig)
{
    loff_t new_pos = 0;
    switch (orig) {
    case 0: /* SEEK_SET: */
        new_pos = offset;
        break;
    case 1: /* SEEK_CUR: */
        new_pos = file->f_pos + offset;
        break;
    case 2: /* SEEK_END: */
        new_pos = MAX_LENGTH - offset;
        break;
    }

    if (new_pos > MAX_LENGTH)
        new_pos = MAX_LENGTH;  // max case
    if (new_pos < 0)
        new_pos = 0;        // min case
    file->f_pos = new_pos;  // This is what we'll use now
    return new_pos;
}

const struct file_operations fib_fops = {
    .owner = THIS_MODULE,
    .read = fib_read,
    .write = fib_write,
    .open = fib_open,
    .release = fib_release,
    .llseek = fib_device_lseek,
};

static int __init init_fib_dev(void)
{
    int rc = 0;

    mutex_init(&fib_mutex);

    // Let's register the device
    // This will dynamically allocate the major number
    rc = alloc_chrdev_region(&fib_dev, 0, 1, DEV_FIBONACCI_NAME);

    if (rc < 0) {
        printk(KERN_ALERT
               "Failed to register the fibonacci char device. rc = %i",
               rc);
        return rc;
    }

    fib_cdev = cdev_alloc();
    if (fib_cdev == NULL) {
        printk(KERN_ALERT "Failed to alloc cdev");
        rc = -1;
        goto failed_cdev;
    }
    fib_cdev->ops = &fib_fops;
    rc = cdev_add(fib_cdev, fib_dev, 1);

    if (rc < 0) {
        printk(KERN_ALERT "Failed to add cdev");
        rc = -2;
        goto failed_cdev;
    }

    fib_class = class_create(THIS_MODULE, DEV_FIBONACCI_NAME);

    if (!fib_class) {
        printk(KERN_ALERT "Failed to create device class");
        rc = -3;
        goto failed_class_create;
    }

    if (!device_create(fib_class, NULL, fib_dev, NULL, DEV_FIBONACCI_NAME)) {
        printk(KERN_ALERT "Failed to create device");
        rc = -4;
        goto failed_device_create;
    }
    return rc;
failed_device_create:
    class_destroy(fib_class);
failed_class_create:
    cdev_del(fib_cdev);
failed_cdev:
    unregister_chrdev_region(fib_dev, 1);
    return rc;
}

static void __exit exit_fib_dev(void)
{
    mutex_destroy(&fib_mutex);
    device_destroy(fib_class, fib_dev);
    class_destroy(fib_class);
    cdev_del(fib_cdev);
    unregister_chrdev_region(fib_dev, 1);
}

module_init(init_fib_dev);
module_exit(exit_fib_dev);
