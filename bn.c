#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bn.h"

#define MAX(x, y) ((x) > (y) ? (x) : (y))
#ifndef SWAP
#define SWAP(x, y)           \
    do {                     \
        typeof(x) __tmp = x; \
        x = y;               \
        y = __tmp;           \
    } while (0)
#endif

#ifndef DIV_ROUNDUP
#define DIV_ROUNDUP(x, len) (((x) + (len) -1) / (len))
#endif


/* count leading zeros of src*/
static int bn_clz(const bn *src)
{
    int cnt = 0;
    for (int i = src->size - 1; i >= 0; i--) {
        if (src->number[i]) {
            // prevent undefined behavior when src = 0
            cnt += __builtin_clz(src->number[i]);
            return cnt;
        } else {
            cnt += 32;
        }
    }
    return cnt;
}

/* count the digits of most significant bit */
static int bn_msb(const bn *src)
{
    return src->size * 32 - bn_clz(src);
}

/*
 * output bn to decimal string
 * Note: the returned string should be freed with the kfree()
 */
char *bn_to_string(const bn *src)
{
    // log10(x) = log2(x) / log2(10) ~= log2(x) / 3.322
    size_t len = (8 * sizeof(int) * src->size) / 3 + 2 + src->sign;
    char *s = malloc(len);
    char *p = s;

    memset(s, '0', len - 1);
    s[len - 1] = '\0';

    /* src.number[0] contains least significant bits */
    for (int i = src->size - 1; i >= 0; i--) {
        /* walk through every bit of bn */
        for (unsigned long long d = 1U << 31; d; d >>= 1) {
            /* binary -> decimal string */
            int carry = !!(d & src->number[i]);
            for (int j = len - 2; j >= 0; j--) {
                s[j] += s[j] - '0' + carry;
                carry = (s[j] > '9');
                if (carry)
                    s[j] -= 10;
            }
        }
    }
    // skip leading zero
    while (p[0] == '0' && p[1] != '\0') {
        p++;
    }
    if (src->sign)
        *(--p) = '-';
    memmove(s, p, strlen(p) + 1);
    return s;
}

/*
 * alloc a bn structure with the given size
 * the value is initialized to +0
 */
bn *bn_alloc(size_t size)
{
    bn *new = (bn *) malloc(sizeof(bn));
    new->number = malloc(sizeof(int) * size);
    new->capacity = size > INIT_ALLOC_SIZE ? size : INIT_ALLOC_SIZE;
    memset(new->number, 0, sizeof(int) * size);
    new->size = size;
    new->sign = 0;
    return new;
}

/*
 * free entire bn data structure
 * return 0 on success, -1 on error
 */
int bn_free(bn *src)
{
    if (src == NULL)
        return -1;
    free(src->number);
    free(src);
    return 0;
}

/*
 * resize bn
 * return 0 on success, -1 on error
 * data lose IS neglected when shinking the size
 */
static int bn_resize(bn *src, size_t size)
{
    if (!src)
        return -1;
    if (size == src->size)
        return 0;
    if (size == 0)  // prevent krealloc(0) = kfree, which will cause problem
        return bn_free(src);
    if (size > src->capacity) {
        src->capacity =
            (size + (ALLOC_CHUNK_SIZE - 1)) & ~(ALLOC_CHUNK_SIZE - 1);
        src->number = realloc(src->number, sizeof(int) * size);
    }
    if (!src->number) {  // realloc fails
        return -1;
    }
    if (size > src->size)
        memset(src->number + src->size, 0, sizeof(int) * (size - src->size));
    src->size = size;
    return 0;
}

/*
 * copy the value from src to dest
 * return 0 on success, -1 on error
 */
int bn_cpy(bn *dest, bn *src)
{
    if (bn_resize(dest, src->size) < 0)
        return -1;
    dest->sign = src->sign;
    memcpy(dest->number, src->number, src->size * sizeof(int));
    return 0;
}

/* swap bn ptr */
void bn_swap(bn *a, bn *b)
{
    bn tmp = *a;
    *a = *b;
    *b = tmp;
}

/* left bit shift on bn (maximun shift 31)
void bn_lshift(bn *src, size_t shift, bn *dest)
{
    size_t z = bn_clz(src);
    shift %= 32;  // only handle shift within 32 bits atm
    if (!shift)
        return;
    if (shift > z) {
        bn_resize(dest, src->size + 1);
    } else {
        bn_resize(dest, src->size);
    }
    for (int i = src->size - 1; i > 0; i--)
        dest->number[i] =
            src->number[i] << shift | src->number[i - 1] >> (32 - shift);
    dest->number[0] = src->number[0] << shift;
}
 right bit shift on bn (maximun shift 31)
void bn_rshift(bn *src, size_t shift)
{
    size_t z = 32 - bn_clz(src);
    shift %= 32;  // only handle shift within 32 bits atm
    if (!shift)
        return;
    for (int i = 0; i < (src->size - 1); i++)
        src->number[i] = src->number[i] >> shift | src->number[i + 1]
                                                       << (32 - shift);
    src->number[src->size - 1] >>= shift;
    if (shift >= z && src->size > 1)
        bn_resize(src, src->size - 1);
}*/

/*
 * compare length
 * return 1 if |a| > |b|
 * return -1 if |a| < |b|
 * return 0 if |a| = |b|
 */
int bn_cmp(const bn *a, const bn *b)
{
    if (a->size > b->size) {
        return 1;
    } else if (a->size < b->size) {
        return -1;
    } else {
        for (int i = a->size - 1; i >= 0; i--) {
            if (a->number[i] > b->number[i])
                return 1;
            if (a->number[i] < b->number[i])
                return -1;
        }
        return 0;
    }
}

/* |c| = |a| + |b| */
static void bn_do_add(const bn *a, const bn *b, bn *c)
{
    // max digits = max(sizeof(a) + sizeof(b)) + 1
    // int d = MAX(bn_msb(a), bn_msb(b));
    int d = MAX(bn_msb(a), bn_msb(b)) + 1;
    d = DIV_ROUNDUP(d, 32) + !d;
    bn_resize(c, d);  // round up, min size = 1

    // unsigned int carry = 0;
    unsigned long long int carry = 0;
    for (int i = 0; i < c->size; i++) {
        unsigned int tmp1 = (i < a->size) ? a->number[i] : 0;
        unsigned int tmp2 = (i < b->size) ? b->number[i] : 0;
        carry += (unsigned long long int) tmp1 + tmp2;
        c->number[i] = carry;
        carry >>= 32;
    }
    if (!c->number[c->size - 1] && c->size > 1)
        bn_resize(c, c->size - 1);
    /*
    for (int i = 0; i < a->size; ++i) {
        unsigned int tmp1 = a->number[i];
        unsigned int tmp2 = b->number[i];
        carry = (tmp2 += carry) < carry;
        carry += (c->number[i] = tmp1 + tmp2) < tmp1;
    }
    if (b->size > a->size) {
        for (int i = a->size; i < b->size; ++i) {
            unsigned int tmp = b->number[i];
            carry = (tmp += carry) < carry;
            c->number[i] = tmp;
        }
    }
    if (carry) {
        c->number[b->size] = carry;
        ++(c->size);
    }*/
}

/*
 * |c| = |a| - |b|
 * Note: |a| > |b| must be true
 */
static void bn_do_sub(const bn *a, const bn *b, bn *c)
{
    // max digits = max(sizeof(a) + sizeof(b))
    int d = MAX(a->size, b->size);
    bn_resize(c, d);

    long long int carry = 0;
    for (int i = 0; i < c->size; i++) {
        unsigned int tmp1 = (i < a->size) ? a->number[i] : 0;
        unsigned int tmp2 = (i < b->size) ? b->number[i] : 0;

        carry = (long long int) tmp1 - tmp2 - carry;
        if (carry < 0) {
            c->number[i] = carry + (1LL << 32);
            carry = 1;
        } else {
            c->number[i] = carry;
            carry = 0;
        }
    }

    d = bn_clz(c) / 32;
    if (d == c->size)
        --d;
    bn_resize(c, c->size - d);
}

/*
 * c = a + b
 * Note: work for c == a or c == b
 */
void bn_add(const bn *a, const bn *b, bn *c)
{
    if (a->sign == b->sign) {  // both positive or negative
        bn_do_add(a, b, c);
        c->sign = a->sign;
    } else {          // different sign
        if (a->sign)  // let a > 0, b < 0
            SWAP(a, b);
        int cmp = bn_cmp(a, b);
        if (cmp > 0) {
            /* |a| > |b| and b < 0, hence c = a - |b| */
            bn_do_sub(a, b, c);
            c->sign = 0;
        } else if (cmp < 0) {
            /* |a| < |b| and b < 0, hence c = -(|b| - |a|) */
            bn_do_sub(b, a, c);
            c->sign = 1;
        } else {
            /* |a| == |b| */
            bn_resize(c, 1);
            c->number[0] = 0;
            c->sign = 0;
        }
    }
}

/*
 * c = a - b
 * Note: work for c == a or c == b
 */
/*
void bn_sub(const bn *a, const bn *b, bn *c)
{*/
/* xor the sign bit of b and let bn_add handle it */
/* bn tmp = *b;
 tmp.sign ^= 1;  // a - b = a + (-b)
 bn_add(a, &tmp, c);
}*/

/* c += x, starting from offset */
/*
static void bn_mult_add(bn *c, int offset, unsigned long long int x)
{
    unsigned long long int carry = 0;
    for (int i = offset; i < c->size; i++) {
        carry += c->number[i] + (x & 0xFFFFFFFF);
        c->number[i] = carry;
        carry >>= 32;
        x >>= 32;
        if (!x && !carry)  // done
            return;
    }
}*/

static unsigned int bn_mult_partial(const unsigned int *num,
                                    unsigned int sz,
                                    unsigned int b,
                                    unsigned int *c)
{
    if (b == 0) {
        return 0;
    }
    unsigned long long carry = 0;
    for (int i = 0; i < sz; ++i) {
        unsigned int high, low;
        unsigned long long product = (unsigned long long) num[i] * b;
        low = product;
        high = product >> 32;
        carry = high + ((low += carry) < carry);
        carry += ((c[i] += low) < low);
    }
    return carry;
} /*
 void bn_sqrt_base(const unsigned int *a, unsigned int sz, unsigned int *c)
 {
     unsigned int *cp = c + 1;
     unsigned int asize = sz - 1;
     const unsigned int *ap = a;
     for (int i = 0; i < asize; ++i) {
         cp[asize - i] = bn_mult_partial(&ap[i+1], asize-i, ap[i], cp);
         cp += 2;
     }
     for (int i = 2 * sz - 1; i > 0; i--)
         c[i] = c[i] << 1 | c[i - 1] >> (31);
     c[0] <<= 1;
     cp = c;
     ap = a;
     asize = sz;
     unsigned long long carry = 0;
     for (int i = 0; i < asize; i++) {
         unsigned int high, low;
     unsigned long long product = (unsigned long long)ap[i] * ap[i];
     low = product;
     high = product >> 32;
         high += (low += carry) < carry;
         high += (cp[0] += low) < low;
         carry = (cp[1] += high) < high;
         cp += 2;
     }
 }*/
void bn_mult_sup(const unsigned int *a,
                 const unsigned int *b,
                 unsigned int a_sz,
                 unsigned int b_sz,
                 unsigned int *c)
{
    if (!a_sz || !b_sz) {
        return;
    }
    for (int i = 0; i < b_sz; ++i) {
        c[a_sz + i] = bn_mult_partial(a, a_sz, b[j], c + i);
    }
}
static unsigned int bn_add_sup(const unsigned *a,
                               const unsigned int *b,
                               unsigned int sz,
                               unsigned int *c)
{
    unsigned int carry = 0;
    for (int i = 0; i < sz; ++i) {
        unsigned int tmp1 = a[i];
        unsigned int tmp2 = b[i];
        carry = (tmp1 += carry) < carry;
        carry += (c[i] = tmp1 + tmp2) < tmp2;
    }

    return carry;
}
unsigned int bn_sub_sup(const unsigned int *a,
                        const unsigned int *b,
                        unsigned int sz,
                        unsigned int *c)
{
    unsigned int borrow = 0;
    for (int i = 0; i < sz; ++i) {
        unsigned int t1 = a[i];
        unsigned int t2 = b[i];
        borrow = (t2 += borrow) < borrow;
        borrow += (c[i] = t1 - t2) > t1;
    }
    return borrow;
}

// Assume the length of two numbers are same. |a| = |b|
void bn_karatsuba(const unsigned int *a,
                  const unsigned int *b,
                  unsigned int sz,
                  unsigned int *c)
{
    unsigned int even_sz = (sz >> 1) << 1;
    unsigned int half_pos = even_sz / 2;
    const unsigned int *x0 = a;
    const unsigned int *x1 = a + half_pos;
    const unsigned int *y0 = b;
    const unsigned int *y1 = b + half_pos;
    unsigned int *cp = c;
    unsigned int *cph = c + even_sz;
    // z0 = x0y0 , z1 = x0y1 + x1y0 , z2 = x1y1
    // z1 = (x0 + x1)(y0 + y1) - z0 - z2 , only one multplication.
    bn_mult_sup(x0, half_pos, y0, half_pos, cp);
    bn_mult_sup(x1, half_pos, y1, half_pos, cph);

    // Calculate the z1 = (x0+x1)(y0+y1) - z0 - z2

    unsigned int *copy_cp =
        (unsigned int *) calloc(even_sz, sizeof(unsigned int));
    for (int i = 0; i < even_sz; ++i) {
        copy_cp[i] = cp[i];
    }
    // Now,add the value from c[half_pos] to c[half_pos + even_sz - 1]
    unsigned int carry = 0;
    carry = bn_add_sup(c + half_pos, cph, even_sz, c + half_pos);
    carry += bn_add_sup(c + half_pos, copy_cp, even_sz, c + half_pos);

    // Calculate the value (x1-x0)
    unsigned int *a_tmp = copy_cp;
    bool prod_neg = bn_cmp(x1, x0);
    if (prod_neg) {
        bn_sub_sup(x0, x1, half_pos, a_tmp);
    } else {
        bn_sub_sup(x1, x0, half_pos, a_tmp);
    }
    // Calculate the value (y0-y1)
    unsigned int *b_tmp = copy_cp + half_pos;
    bool prod_neg1 = bn_cmp(x1, x0);
    if (prod_neg1) {
        bn_sub_sup(y1, y0, half_pos, b_tmp);
    } else {
        bn_sub_sup(y0, y1, half_pos, b_tmp);
    }

    // Calculate the value |(x1-x0)|*|(y0-y1)|
    tmp = (unsigned int *) (calloc(even_sz, sizeof(unsigned int)));
    if (half_pos > KARATSUBA_MUL_THRESHOLD) {
        bn_karatsuba(a_tmp, b_tmp, half_pos, tmp);
    } else {
        bn_mult_sup(a_tmp, b_tmp, half_pos, half_pos, tmp);
    }
    free(a_tmp);

    // Add the (x1-x0)*(y0-y1) to the c[half_pos ... half_pos + even_sz -1]
    // Also , considet the sign of the result (x1-x0)*(y0-y1)
    if (prod_neg ^ prod_neg1)
        carry -= bn_sub_sup(c + half_pos, tmp, even_sz, c + half_pos);
    else
        carry += bn_add_sup(c + half_pos, tmp, even_sz, c + half_pos);

    free(tmp);

    // Add the carry to the c[half_pos+even_sz ... 2*even_sz - 1]
    for (int i = half_pos + even_size; i < 2 * even_sz; ++i) {
        unsigned int temp = c[i];
        carry = (temp += carry) < carry;
        c[i] = carry;
    }

    // If sz is the odd
    if (sz & 1) {
        c[even_sz * 2] = bn_mult_partial();
        c[even_sz * 2 + 1] = bn_mult_partial();
    }
}

/*
 * c = a x b
 * Note: work for c == a or c == b
 * using the simple quadratic-time algorithm (long multiplication)
 */
void bn_mult(const bn *a, const bn *b, bn *c)
{
    // max digits = sizeof(a) + sizeof(b))
    int d = bn_msb(a) + bn_msb(b);
    d = DIV_ROUNDUP(d, 32) + !d;  // round up, min size = 1
    bn *tmp;
    /* make it work properly when c == a or c == b */
    if (c == a || c == b) {
        tmp = c;  // save c
        c = bn_alloc(d);
    } else {
        tmp = NULL;
        for (int i = 0; i < c->size; i++)
            c->number[i] = 0;
        bn_resize(c, d);
    }
    /*
    for (int i = 0; i < a->size; i++) {
        for (int j = 0; j < b->size; j++) {
            unsigned long long int carry = 0;
            carry = (unsigned long long int) a->number[i] * b->number[j];
            bn_mult_add(c, i + j, carry);
        }
    }
    */

    for (int i = 0; i < a->size; ++i) {
        c->number[a->size + i] =
            bn_mult_partial(a->number, a->size, b->number[i], c->number + i);
    }
    c->sign = a->sign ^ b->sign;

    if (tmp) {
        bn_cpy(tmp, c);  // restore c
        bn_free(c);
    }
}

/* calc n-th Fibonacci number and save into dest */
/*void bn_fib(bn *dest, unsigned int n)
{
    bn_resize(dest, 1);
    if (n < 2) {  // Fib(0) = 0, Fib(1) = 1
        dest->number[0] = n;
        return;
    }
    bn *a = bn_alloc(1);
    bn *b = bn_alloc(1);
    dest->number[0] = 1;
    for (unsigned int i = 1; i < n; i++) {
        bn_cpy(b, dest);
        bn_add(dest, a, dest);
        bn_swap(a, b);
    }
    bn_free(a);
    bn_free(b);
}*/
